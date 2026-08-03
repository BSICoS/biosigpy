"""Slope-range ECG-derived respiration."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def sloperange(
    decg: ArrayLike,
    r_wave_times: ArrayLike,
    sampling_frequency: float,
) -> np.ndarray:
    """Estimate respiration from derivative ECG morphology around R waves.

    Parameters
    ----------
    decg : array_like
        One-dimensional derivative ECG signal with at least two finite real
        samples.
    r_wave_times : array_like
        Finite R-wave occurrence times in seconds. Values must form a
        non-empty, strictly increasing one-dimensional sequence and map inside
        the derivative ECG sample grid.
    sampling_frequency : float
        Positive finite sampling frequency in Hz.

    Returns
    -------
    numpy.ndarray
        One-dimensional ECG-derived respiration amplitude vector aligned with
        ``r_wave_times``. Beats whose analysis windows cross a signal boundary
        are represented by ``NaN``.

    Raises
    ------
    TypeError
        If an input expected to be numeric contains non-numeric or complex
        data.
    ValueError
        If an input has an invalid shape, length, finite-value constraint,
        order, sample-grid mapping, or sampling frequency.

    Notes
    -----
    This function implements the normative Biosiglib ``ecg.sloperange``
    output. Diagnostic slope arrays exposed by some implementations are not
    part of the public contract and are intentionally not returned.

    Examples
    --------
    >>> import numpy as np
    >>> from biosigpy.ecg import sloperange
    >>> decg = np.zeros(40)
    >>> decg[[9, 19, 29]] = [3.0, 6.0, 2.0]
    >>> decg[[13, 23, 33]] = [-2.0, -1.0, -4.0]
    >>> sloperange(decg, [0.1, 0.2, 0.3], 100.0)
    array([5., 7., 6.])
    """

    derivative_ecg = as_real_vector(decg, name="decg")
    if derivative_ecg.size < 2:
        raise ValueError("decg must contain at least two samples")
    if not np.all(np.isfinite(derivative_ecg)):
        raise ValueError("decg must contain only finite values")

    times = as_real_vector(r_wave_times, name="r_wave_times")
    if times.size == 0:
        raise ValueError("r_wave_times must not be empty")
    if not np.all(np.isfinite(times)):
        raise ValueError("r_wave_times must contain only finite values")
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("r_wave_times must be strictly increasing")

    fs = as_positive_real_scalar(
        sampling_frequency, name="sampling_frequency"
    )
    with np.errstate(over="ignore", invalid="ignore"):
        sample_positions = times * fs
    if not np.all(np.isfinite(sample_positions)):
        raise ValueError("r_wave_times cannot be mapped to the sample grid")

    rounded_positions = _round_half_away_from_zero(sample_positions)
    if np.any(rounded_positions < 0) or np.any(
        rounded_positions >= derivative_ecg.size
    ):
        raise ValueError("r_wave_times must map inside the decg sample grid")
    r_wave_samples = rounded_positions.astype(np.intp)

    short_window = int(np.floor(fs * 0.015 + 0.5))
    long_window = int(np.floor(fs * 0.05 + 0.5))
    edr = np.full(times.shape, np.nan, dtype=np.float64)

    for beat_index, sample in enumerate(r_wave_samples):
        sample = int(sample)
        upslope_start = sample - long_window + 1
        upslope_stop = sample + short_window + 1
        downslope_start = sample - short_window
        downslope_stop = sample + long_window

        windows_are_empty = (
            upslope_start >= upslope_stop
            or downslope_start >= downslope_stop
        )
        windows_cross_boundary = (
            upslope_start < 0
            or downslope_start < 0
            or upslope_stop > derivative_ecg.size
            or downslope_stop > derivative_ecg.size
        )
        if windows_are_empty or windows_cross_boundary:
            continue

        maximum_upslope = np.max(
            derivative_ecg[upslope_start:upslope_stop]
        )
        minimum_downslope = np.min(
            derivative_ecg[downslope_start:downslope_stop]
        )
        edr[beat_index] = maximum_upslope - minimum_downslope

    return edr


def _round_half_away_from_zero(values: np.ndarray) -> np.ndarray:
    return np.copysign(np.floor(np.abs(values) + 0.5), values)
