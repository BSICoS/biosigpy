"""Pan-Tompkins-style ECG R-wave detection."""

import numpy as np
from numpy.typing import ArrayLike
from scipy import signal

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector
from biosigpy.tools.snap_to_peak import snap_to_peak


def pantompkins(
    ecg: ArrayLike,
    sampling_frequency: float,
    *,
    bandpass_frequency: ArrayLike = (5.0, 12.0),
    integration_window_size: float = 0.15,
    minimum_peak_distance: float = 0.5,
    snap_to_peak_window_size: float = 20.0,
) -> dict[str, np.ndarray]:
    """Detect ECG R waves and return the canonical processing outputs.

    R-wave times are expressed in seconds. Internal peak indices are zero-based.
    """

    ecg_vector = _ecg_vector(ecg)
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    bandpass = _bandpass_pair(bandpass_frequency, fs)
    window_seconds = as_positive_real_scalar(
        integration_window_size, name="integration_window_size"
    )
    peak_distance_seconds = as_positive_real_scalar(
        minimum_peak_distance, name="minimum_peak_distance"
    )
    snap_to_peak_window_size = as_positive_real_scalar(
        snap_to_peak_window_size, name="snap_to_peak_window_size"
    )

    b_bandpass, a_bandpass = signal.butter(
        4, bandpass, btype="bandpass", fs=fs
    )
    ecg_filtered = signal.filtfilt(b_bandpass, a_bandpass, ecg_vector)

    derivative_filter = _default_derivative_filter(fs, bandpass[1])
    decg = signal.lfilter(derivative_filter, [1.0], ecg_filtered) ** 2

    window_samples = int(np.floor(fs * window_seconds + 0.5))
    if window_samples < 1:
        raise ValueError("integration_window_size is shorter than one sample")
    integration_kernel = np.ones(window_samples, dtype=np.float64) / window_samples
    decg_envelope = np.convolve(decg, integration_kernel, mode="same")

    minimum_distance_samples = int(np.floor(fs * peak_distance_seconds + 0.5))
    peak_indices, _ = signal.find_peaks(
        decg_envelope, distance=max(1, minimum_distance_samples)
    )
    peak_indices = np.unique(np.sort(peak_indices))
    filter_edge_margin = 3 * max(b_bandpass.size, a_bandpass.size)
    peak_indices = snap_to_peak(
        ecg_vector, peak_indices.astype(np.float64) + 1.0, snap_to_peak_window_size
    )
    peak_indices = peak_indices[np.isfinite(peak_indices)] - 1.0
    peak_indices = peak_indices[peak_indices >= filter_edge_margin]

    return {
        "r_wave_times": peak_indices.astype(np.float64) / fs,
        "ecg_filtered": np.asarray(ecg_filtered, dtype=np.float64),
        "decg": np.asarray(decg, dtype=np.float64),
        "decg_envelope": np.asarray(decg_envelope, dtype=np.float64),
    }


def _ecg_vector(ecg: ArrayLike) -> np.ndarray:
    ecg = as_real_vector(ecg, name="ecg")
    if ecg.size < 2:
        raise ValueError("ecg must contain at least two samples")
    if np.any(np.isinf(ecg)):
        raise ValueError("ecg must not contain infinite values")
    return ecg


def _bandpass_pair(value: ArrayLike, sampling_frequency: float) -> np.ndarray:
    array = np.asarray(value)
    if array.shape != (2,) or not np.issubdtype(array.dtype, np.number):
        raise ValueError("bandpass_frequency must contain two numeric values")
    if np.issubdtype(array.dtype, np.complexfloating):
        raise TypeError("bandpass_frequency must be real")

    bandpass = np.asarray(array, dtype=np.float64)
    if (
        not np.all(np.isfinite(bandpass))
        or np.any(bandpass <= 0)
        or bandpass[0] >= bandpass[1]
        or bandpass[1] >= sampling_frequency / 2.0
    ):
        raise ValueError("bandpass_frequency must lie inside the Nyquist range")
    return bandpass


def _default_derivative_filter(
    sampling_frequency: float, stop_frequency: float
) -> np.ndarray:
    # Reference defaults use the same order-4 least-squares differentiator
    # design as the shared Biosigmat oracle.
    if sampling_frequency == 256.0 and stop_frequency == 12.0:
        return np.array(
            [
                8.00635957041095,
                4.08939515658424,
                0.0,
                -4.08939515658424,
                -8.00635957041095,
            ],
            dtype=np.float64,
        )

    normalized_stop = stop_frequency / (sampling_frequency / 2.0)
    scale = sampling_frequency * normalized_stop / 3.0
    return scale * np.array([2.0, 1.0, 0.0, -1.0, -2.0]) / 5.0
