"""Detect artifact intervals from windowed Hjorth descriptors."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools import hjorth
from biosigpy.tools._validation import as_integer_scalar, as_positive_real_scalar, as_real_vector


def hjorth_artifacts(
    signal: ArrayLike, sampling_frequency: float, window_duration: float,
    step_duration: float, margins: ArrayLike, *, min_segment_separation: float = 0,
    median_filter_order: int = 300, negative: bool = False, plotflag: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a binary artifact mask and inclusive time intervals in seconds.

    Parameters
    ----------
    signal : array_like
        Nonempty real signal with finite samples or NaN.
    sampling_frequency : float
        Finite positive sampling frequency in Hz.
    window_duration, step_duration : float
        Positive durations in seconds, converted to samples by flooring.
        Only complete windows are analyzed; at least three samples per window
        and one sample per step are required.
    margins : array_like, shape (3, 2)
        Nonnegative absolute lower and upper margins around the median
        baselines for activity, mobility and complexity, respectively.
    min_segment_separation : float, optional
        Merge uncovered gaps strictly shorter than this duration in seconds.
        The default zero merges only overlapping or contiguous intervals.
    median_filter_order : int, optional
        Median window length in analysis windows, default 300. Boundaries
        truncate the median window and NaNs are omitted.
    negative : bool, optional
        Invert the threshold decision. Windows containing NaN remain marked.
    plotflag : bool, optional
        Plot the signal, detected artifacts, Hjorth parameters and thresholds
        in four panels, as in MATLAB. Requires Matplotlib. Default is False.

    Returns
    -------
    artifact_mask : ndarray
        One-dimensional boolean mask aligned with signal. Samples outside
        full analysis windows are unmarked unless included in a merged gap.
    artifact_intervals : ndarray, shape (n, 2)
        Times of the first and last marked sample in each merged interval,
        with the first input sample at time zero. No artifacts yields (0, 2).
    """
    # Parse and validate inputs
    signal = as_real_vector(signal, name="signal")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    seg = as_positive_real_scalar(window_duration, name="window_duration")
    step = as_positive_real_scalar(step_duration, name="step_duration")
    margins = np.asarray(margins)
    if margins.dtype.kind not in "iuf" or margins.shape != (3, 2):
        raise ValueError("margins must be a real 3-by-2 matrix")
    if not np.all(np.isfinite(margins)) or np.any(margins < 0):
        raise ValueError("margins must be finite and nonnegative")
    separation = np.asarray(min_segment_separation)
    if separation.ndim != 0 or separation.dtype.kind not in "iuf":
        raise ValueError("min_segment_separation must be a real scalar")
    min_segment_separation = float(separation)
    if not np.isfinite(min_segment_separation) or min_segment_separation < 0:
        raise ValueError("min_segment_separation must be finite and nonnegative")
    median_filter_order = as_integer_scalar(median_filter_order, name="median_filter_order")
    if median_filter_order < 1 or not isinstance(negative, (bool, np.bool_)):
        raise ValueError("median_filter_order must be positive and negative must be boolean")
    if not isinstance(plotflag, (bool, np.bool_)):
        raise ValueError("plotflag must be boolean")
    if signal.size == 0 or np.any(np.isinf(signal)):
        raise ValueError("signal must be nonempty and contain finite samples or NaN")

    # Extract margin values
    margin_h0_low, margin_h0_up = margins[0]
    margin_h1_low, margin_h1_up = margins[1]
    margin_h2_low, margin_h2_up = margins[2]

    # Initialize Hjorth parameters
    n_window = int(np.floor(seg * fs))
    n_step = int(np.floor(step * fs))
    if n_window < 3 or n_step < 1:
        raise ValueError("require at least three samples per window and one per step")
    if signal.size < n_window:
        raise ValueError("signal must contain a complete window")
    n_segments = (signal.size - n_window) // n_step
    h0 = np.zeros(n_segments + 1)
    h1 = np.zeros(n_segments + 1)
    h2 = np.zeros(n_segments + 1)
    invalid_windows = np.zeros(n_segments + 1, dtype=bool)

    # Subtract the global mean, retaining missing samples
    valid_samples = signal[~np.isnan(signal)]
    sfilt = signal - (np.mean(valid_samples) if valid_samples.size else np.nan)

    # Compute Hjorth parameters
    for k in range(n_segments + 1):
        indexes = slice(k * n_step, k * n_step + n_window)
        invalid_windows[k] = np.any(np.isnan(signal[indexes]))
        h0[k], h1[k], h2[k] = hjorth(sfilt[indexes], fs)

    # Compute thresholds
    baseline_h0 = _median_baseline(h0, median_filter_order)
    baseline_h1 = _median_baseline(h1, median_filter_order)
    baseline_h2 = _median_baseline(h2, median_filter_order)
    threshold_h0_low = baseline_h0 - margin_h0_low
    threshold_h0_up = baseline_h0 + margin_h0_up
    threshold_h1_low = baseline_h1 - margin_h1_low
    threshold_h1_up = baseline_h1 + margin_h1_up
    threshold_h2_low = baseline_h2 - margin_h2_low
    threshold_h2_up = baseline_h2 + margin_h2_up
    threshold_h0_low[threshold_h0_low <= 0] = 0.0001

    # Look for artifact segments
    is_artifact = ((h2 > threshold_h2_up) | (h2 < threshold_h2_low) |
                   (h1 > threshold_h1_up) | (h1 < threshold_h1_low) |
                   (h0 > threshold_h0_up) | (h0 < threshold_h0_low))
    if negative:
        is_artifact = ~is_artifact
    detected = np.flatnonzero(is_artifact | invalid_windows)

    artifact_mask = np.zeros(signal.size, dtype=bool)
    artifact_intervals = np.empty((0, 2))
    if detected.size:
        # Convert detected windows to sample intervals and merge in signal order
        artifact_segments = np.column_stack((detected * n_step, detected * n_step + n_window - 1))
        artifact_intervals = [artifact_segments[0].tolist()]
        for start, end in artifact_segments[1:]:
            gap_samples = start - artifact_intervals[-1][1] - 1
            if gap_samples <= 0 or gap_samples / fs < min_segment_separation:
                artifact_intervals[-1][1] = max(artifact_intervals[-1][1], int(end))
            else:
                artifact_intervals.append([int(start), int(end)])
        for start, end in artifact_intervals:
            artifact_mask[start:end + 1] = True
        artifact_intervals = np.asarray(artifact_intervals, dtype=float) / fs

    if plotflag:
        _plot_results(signal, fs, h0, h1, h2, threshold_h0_up, threshold_h0_low,
                      threshold_h1_up, threshold_h1_low, threshold_h2_up,
                      threshold_h2_low, artifact_intervals)

    return artifact_mask, artifact_intervals


def _median_baseline(values: np.ndarray, order: int) -> np.ndarray:
    """Match medfilt1 with truncated boundaries and omitted NaNs."""
    baseline = np.full(values.size, np.nan)
    before = order // 2
    after = (order - 1) // 2
    for k in range(values.size):
        window = values[max(0, k-before):min(values.size, k+after+1)]
        window = window[~np.isnan(window)]
        if window.size:
            baseline[k] = np.median(window)
    return baseline



def _plot_results(signal, fs, h0, h1, h2, threshold_h0_up, threshold_h0_low,
                  threshold_h1_up, threshold_h1_low, threshold_h2_up,
                  threshold_h2_low, detections):
    """Plot the same signal and descriptor panels as MATLAB plotResults."""
    import matplotlib.pyplot as plt

    time = np.linspace(0, (signal.size - 1) / fs, h1.size)
    signal_time = np.arange(signal.size) / fs
    figure, axes = plt.subplots(4, 1, sharex=True, figsize=(12, 9))

    axes[0].plot(signal_time, signal)
    for start, end in detections:
        inside = (signal_time >= start) & (signal_time <= end)
        axes[0].plot(signal_time[inside], signal[inside], "r")
    axes[0].set_ylabel("Input")

    parameters = (h0, h1, h2)
    lower_thresholds = (threshold_h0_low, threshold_h1_low, threshold_h2_low)
    upper_thresholds = (threshold_h0_up, threshold_h1_up, threshold_h2_up)
    labels = ("Activity", "Mobility", "Complexity")
    for axis, values, lower, upper, label in zip(
        axes[1:], parameters, lower_thresholds, upper_thresholds, labels
    ):
        axis.plot(time, values)
        axis.plot(time, upper, "k", linewidth=2)
        axis.plot(time, lower, "k", linewidth=2)
        axis.set_ylabel(label)
        finite_lower = lower[np.isfinite(lower)]
        finite_upper = upper[np.isfinite(upper)]
        if finite_lower.size and finite_upper.size:
            low, high = np.min(finite_lower), np.max(finite_upper)
            if low < high:
                axis.set_ylim(low, high)
    for axis in axes:
        axis.grid(True)
    if time[0] < time[-1]:
        axes[-1].set_xlim(time[0], time[-1])
    axes[-1].set_xlabel("t [sec]")
    figure.tight_layout()
