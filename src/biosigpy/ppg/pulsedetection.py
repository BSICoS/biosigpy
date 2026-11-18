"""Detect maximum upslopes in a filtered PPG derivative."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools import localmax
from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def pulsedetection(
    dppg: ArrayLike,
    sampling_frequency: float,
    *,
    method: str = "adaptive",
    adaptive_alpha_amp: float = 0.2,
    adaptive_refract_period: float = 0.15,
    adaptive_tau_rr: float = 1.0,
    reset_after: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return refined pulse maximum-upslope times and an adaptive threshold.

    Parameters
    ----------
    dppg : array_like
        Nonempty real PPG derivative after low-pass differentiation and any
        filter-delay compensation. NaNs are allowed; infinities are invalid.
    sampling_frequency : float
        Finite positive sampling frequency in Hz.
    method : str, optional
        Detection algorithm: ``"adaptive"`` (default), case-insensitive.
    adaptive_alpha_amp : float, optional
        Positive amplitude multiplier defining the threshold floor, default 0.2.
    adaptive_refract_period : float, optional
        Positive refractory threshold-hold duration in seconds, default 0.15.
    adaptive_tau_rr : float, optional
        Positive decay duration after the refractory period, expressed as a
        fraction of the estimated pulse interval. Default is 1.0.
    reset_after : float or None, optional
        Positive seconds without a raw detection before restarting adaptation.
        Missing spans count toward elapsed time; initialization resumes when
        signal is available. None (default) or positive infinity disables resets.

    Returns
    -------
    t_d : ndarray
        One-dimensional array of quadratically refined detection times in
        seconds, with the first input sample at time zero. Undefined refined
        times are omitted. No detections yields an empty array.
    threshold : ndarray
        One-dimensional adaptive threshold aligned with every input sample.

    Notes
    -----
    Adaptive state is preserved throughout the signal, including NaN spans,
    unless the optional inactivity reset is enabled. Missing samples do not
    satisfy crossings, but a crossing region can span them. Raw detections
    update adaptation even if missing neighbors prevent time refinement.

    A finite ``reset_after`` can help after a device gain reduction leaves
    pulses below the previous threshold floor. It is optional and its duration
    depends on the application. This function does not perform preprocessing.
    """
    # Parse and validate inputs
    dppg = as_real_vector(dppg, name="dppg")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    if dppg.size == 0 or np.any(np.isinf(dppg)):
        raise ValueError("dppg must be nonempty and contain finite samples or NaN")
    if not isinstance(method, str) or method.lower() != "adaptive":
        raise ValueError("method must be 'adaptive'")
    alpha_amp = as_positive_real_scalar(adaptive_alpha_amp, name="adaptive_alpha_amp")
    refract_period = as_positive_real_scalar(
        adaptive_refract_period, name="adaptive_refract_period"
    )
    tau_rr = as_positive_real_scalar(adaptive_tau_rr, name="adaptive_tau_rr")
    if reset_after is None:
        reset_after = np.inf
    else:
        reset_value = np.asarray(reset_after)
        if not (reset_value.ndim == 0 and reset_value.dtype.kind in "iuf"
                and np.isposinf(reset_value)):
            reset_after = as_positive_real_scalar(reset_after, name="reset_after")
        else:
            reset_after = np.inf

    # Detect pulses while preserving adaptive state across search blocks
    detection_samples, threshold = _adaptive_threshold(
        dppg, fs, alpha_amp, refract_period, tau_rr, reset_after
    )

    # Remove duplicates and refine detection times
    detection_samples = np.unique(detection_samples)
    t_d = np.empty(0)
    if detection_samples.size:
        left = dppg[detection_samples - 1]
        center = dppg[detection_samples]
        right = dppg[detection_samples + 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            delta = (left - right) / (2 * (left - 2 * center + right))
        t_d = (detection_samples + delta) / fs

    # Remove NaN values from refined detection times
    t_d = t_d[~np.isnan(t_d)]
    return t_d, threshold


@dataclass
class _ThresholdProfile:
    decay_start: int
    amplitude: float
    slope: float
    floor_start: int
    floor_value: float


def _adaptive_threshold(dppg, fs, alpha_amp, refract_period, tau_rr, reset_after):
    # Algorithm constants
    initial_hr_estimate = 80  # bpm
    min_peaks_for_rr_estimation = 4
    n_rr_estimation = 3
    n_amplitude_estimation = 3
    search_block_length = 4096

    # Prepare input signal
    signal_length = dppg.size
    refract_period_samples = _round_samples(refract_period * fs)
    reset_samples = int(np.ceil(reset_after * fs)) if np.isfinite(reset_after) else np.inf

    # Initialize threshold using initial signal window
    initial_rr_samples = _round_samples(60 / initial_hr_estimate * fs)
    profile, first_signal_index = _initialize_threshold(
        dppg, fs, alpha_amp, initial_rr_samples, 0, search_block_length
    )
    threshold = np.full(signal_length, np.nan)
    threshold_write_start = 0
    current_rr_samples = initial_rr_samples
    reset_index = np.inf
    if first_signal_index is not None:
        reset_index = first_signal_index + reset_samples

    # Initialize pulse detection variables
    detection_samples = []
    history_start = 0

    # Main pulse detection loop
    current_index = 0
    while True:
        search_end = int(min(signal_length - 1, reset_index - 1))

        # Find next threshold crossing (signal rises above threshold)
        up_crossing = _find_crossing(
            dppg, profile, current_index, search_end, search_block_length, True
        )

        # Find where signal drops back below threshold
        down_crossing = None
        if up_crossing is not None:
            down_crossing = _find_crossing(
                dppg, profile, up_crossing, search_end, search_block_length, False
            )

        if up_crossing is None or down_crossing is None:
            if reset_index >= signal_length:
                # No more complete upward crossings - end detection
                break

            # Reinitialize only when the caller has enabled the inactivity limit.
            next_signal_index = _find_available_sample(
                dppg, int(reset_index), search_block_length
            )
            if next_signal_index is None:
                break
            threshold[threshold_write_start:next_signal_index] = _threshold_values(
                profile, np.arange(threshold_write_start, next_signal_index)
            )
            profile, first_signal_index = _initialize_threshold(
                dppg, fs, alpha_amp, initial_rr_samples, next_signal_index,
                search_block_length,
            )
            threshold_write_start = next_signal_index
            current_rr_samples = initial_rr_samples
            history_start = len(detection_samples)
            current_index = next_signal_index
            reset_index = first_signal_index + reset_samples
            continue

        # Find peak within the crossing region, including its preceding neighbor
        peak_search_start = max(0, up_crossing - 1)
        _, relative_peak_index = localmax(dppg[peak_search_start:down_crossing + 1])
        if np.isnan(relative_peak_index):
            # No peak found - move to next index
            current_index = down_crossing + 1
            continue

        # Find peak index
        peak_index = peak_search_start + int(relative_peak_index)

        # Calculate peak amplitude using adaptive estimation
        num_history_pulses = len(detection_samples) - history_start
        current_peak = np.nanmax(dppg[up_crossing:down_crossing + 1])
        if num_history_pulses <= min_peaks_for_rr_estimation:
            current_amplitude = current_peak
        else:
            # Use median of recent peaks for more robust estimation
            recent_peaks = dppg[detection_samples[-(n_amplitude_estimation + 1):]]
            current_amplitude = np.median(np.r_[recent_peaks, current_peak])

        # Store detected pulse
        detection_samples.append(peak_index)
        num_history_pulses = len(detection_samples) - history_start
        reset_index = peak_index + reset_samples

        # Update RR interval estimate based on detected pulses
        if num_history_pulses >= n_rr_estimation + 1:
            # Use median of recent RR intervals for robust estimation
            recent_rr_intervals = np.diff(detection_samples[-(n_rr_estimation + 1):])
            current_rr_samples = _round_samples(np.median(recent_rr_intervals))
        elif num_history_pulses >= 2:
            # Use mean when insufficient data for robust median
            current_rr_samples = _round_samples(
                np.mean(np.diff(detection_samples[history_start:]))
            )

        # Set refractory period after pulse detection
        refractory_end = min(peak_index + refract_period_samples, signal_length - 1)
        current_index = refractory_end

        # Finalize the preceding threshold without filling the entire future signal
        threshold[threshold_write_start:peak_index] = _threshold_values(
            profile, np.arange(threshold_write_start, peak_index)
        )

        # Set threshold to peak amplitude during refractory period
        threshold[peak_index:refractory_end + 1] = current_amplitude

        # Calculate threshold fall-off parameters
        fall_threshold = current_amplitude * alpha_amp

        # Apply amplitude-based correction if we have sufficient pulse history
        if num_history_pulses >= n_amplitude_estimation + 1:
            recent_amplitudes = dppg[detection_samples[-(n_amplitude_estimation + 1):-1]]
            estimated_amplitude = np.median(recent_amplitudes)

            # Correct for amplitude outliers (peaks that are too high)
            if current_amplitude >= 2 * estimated_amplitude:
                fall_threshold = alpha_amp * estimated_amplitude
                current_amplitude = estimated_amplitude

        # Calculate threshold decay profile after refractory period
        fall_duration = _round_samples(tau_rr * current_rr_samples)
        with np.errstate(divide="ignore", invalid="ignore"):
            slope = np.divide(current_amplitude - fall_threshold, fall_duration)
        profile = _ThresholdProfile(
            refractory_end, current_amplitude, slope,
            refractory_end + fall_duration + 1, fall_threshold,
        )
        threshold_write_start = refractory_end

    # Complete the final threshold profile
    threshold[threshold_write_start:] = _threshold_values(
        profile, np.arange(threshold_write_start, signal_length)
    )
    return np.asarray(detection_samples, dtype=int), threshold


def _initialize_threshold(dppg, fs, alpha_amp, initial_rr_samples, start_index, block_length):
    """Initialize the amplitude estimate and linear decay."""
    initial_window_duration = 10  # seconds
    amplitude_multiplier = 3
    first_signal_index = _find_available_sample(dppg, start_index, block_length)
    initial_threshold = np.nan
    if first_signal_index is not None:
        init_window_end = min(
            first_signal_index + _round_samples(initial_window_duration * fs), dppg.size - 1
        )
        initial_window = dppg[first_signal_index:init_window_end + 1]
        nonnegative_samples = initial_window[initial_window >= 0]
        if nonnegative_samples.size:
            initial_threshold = amplitude_multiplier * np.mean(nonnegative_samples)

    with np.errstate(divide="ignore", invalid="ignore"):
        slope = np.divide(initial_threshold * (1 - alpha_amp), initial_rr_samples)
    floor_start = dppg.size
    if start_index + initial_rr_samples < dppg.size - 1:
        floor_start = start_index + initial_rr_samples
    profile = _ThresholdProfile(
        start_index, initial_threshold, slope, floor_start, alpha_amp * initial_threshold
    )
    return profile, first_signal_index


def _threshold_values(profile, indices):
    """Evaluate only the requested part of the current threshold profile."""
    with np.errstate(invalid="ignore"):
        values = profile.amplitude - profile.slope * (indices - profile.decay_start)
    values[indices >= profile.floor_start] = profile.floor_value
    return values


def _find_crossing(dppg, profile, start_index, end_index, block_length, upward):
    """Search to the end using bounded temporary arrays and no resets."""
    for block_start in range(start_index, end_index + 1, block_length):
        block_end = min(block_start + block_length, end_index + 1)
        values = _threshold_values(profile, np.arange(block_start, block_end))
        signal = dppg[block_start:block_end]
        matches = np.flatnonzero(signal > values if upward else signal < values)
        if matches.size:
            return block_start + int(matches[0])
    return None


def _find_available_sample(dppg, start_index, block_length):
    """Skip missing samples without a mask for the entire remaining signal."""
    for block_start in range(start_index, dppg.size, block_length):
        matches = np.flatnonzero(~np.isnan(dppg[block_start:block_start + block_length]))
        if matches.size:
            return block_start + int(matches[0])
    return None


def _round_samples(value):
    """Round nonnegative sample counts with half-way values upward."""
    return int(np.floor(value + 0.5))
