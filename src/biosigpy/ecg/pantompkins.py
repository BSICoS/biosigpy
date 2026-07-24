"""Pan-Tompkins-style ECG R-wave detection."""

import numpy as np
from numpy.typing import ArrayLike
from scipy import signal

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector
from biosigpy.tools.lpd_filter import lpd_filter
from biosigpy.tools.nan_filter import nan_filter
from biosigpy.tools.nan_filtfilt import nan_filtfilt
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
    """Detect ECG R waves with the Biosiglib Pan-Tompkins workflow.

    Parameters
    ----------
    ecg : array_like
        One-dimensional real numeric ECG signal with at least two samples.
        Infinite values are invalid.
    sampling_frequency : float
        Positive sampling frequency in Hz.
    bandpass_frequency : array_like, optional
        Two-element bandpass frequency vector in Hz. Values must be positive,
        strictly increasing, and inside the Nyquist range.
    integration_window_size : float, optional
        Moving-integration window length in seconds.
    minimum_peak_distance : float, optional
        Minimum distance between envelope peaks in seconds.
    snap_to_peak_window_size : float, optional
        Search radius in samples used when refining detections to local ECG
        maxima.

    Returns
    -------
    dict[str, numpy.ndarray]
        Dictionary with ``r_wave_times`` in seconds, ``ecg_filtered``,
        ``decg_squared``, and ``decg_envelope``.

    Raises
    ------
    TypeError
        If numeric inputs are not real numeric data.
    ValueError
        If inputs are outside the accepted shape, range, or positivity
        constraints.

    Notes
    -----
    This function follows the Biosiglib ``ecg.pantompkins`` specification. It
    returns R-wave times and intermediate processing signals. Bandpass and
    derivative filtering use the NaN-aware public filtering tools with
    ``max_gap=0``. Peak refinement uses
    :func:`biosigpy.tools.snap_to_peak.snap_to_peak`, and derivative filter
    design uses :func:`biosigpy.tools.lpd_filter.lpd_filter`.

    Examples
    --------
    >>> import numpy as np
    >>> from biosigpy.ecg import pantompkins
    >>> fs = 256.0
    >>> time = np.arange(0.0, 2.0, 1.0 / fs)
    >>> ecg = np.sin(2.0 * np.pi * 1.2 * time)
    >>> outputs = pantompkins(ecg, fs)
    >>> "r_wave_times" in outputs
    True
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
    ecg_filtered = nan_filtfilt(
        b_bandpass,
        a_bandpass,
        ecg_vector,
        max_gap=0,
    )

    derivative_filter, _ = lpd_filter(fs, bandpass[1], order=4)
    decg = nan_filter(
        derivative_filter,
        [1.0],
        ecg_filtered,
        max_gap=0,
    )
    decg_squared = decg**2

    window_samples = int(np.floor(fs * window_seconds + 0.5))
    if window_samples < 1:
        raise ValueError("integration_window_size is shorter than one sample")
    integration_kernel = np.ones(window_samples, dtype=np.float64) / window_samples
    decg_envelope = np.convolve(decg_squared, integration_kernel, mode="same")

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
        "decg_squared": np.asarray(decg_squared, dtype=np.float64),
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
