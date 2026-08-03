"""Slope-range ECG-derived respiration."""

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


__all__ = ["SlopeRangeResult", "sloperange"]


class SlopeRangeResult(NamedTuple):
    """Named, unpackable result of :func:`sloperange`.

    Attributes
    ----------
    edr : numpy.ndarray
        ECG-derived respiration amplitudes aligned with R-wave times, in the
        same arbitrary amplitude unit as ``decg``.
    upslopes : numpy.ndarray
        Signal-aligned derivative ECG values inside complete upslope windows,
        with ``NaN`` elsewhere.
    downslopes : numpy.ndarray
        Signal-aligned derivative ECG values inside complete downslope
        windows, with ``NaN`` elsewhere.
    upslope_max_positions : numpy.ndarray
        Zero-based positions of selected upslope maxima aligned with R-wave
        times, in samples, with ``NaN`` for incomplete beats.
    downslope_min_positions : numpy.ndarray
        Zero-based positions of selected downslope minima aligned with R-wave
        times, in samples, with ``NaN`` for incomplete beats.
    """

    edr: np.ndarray
    upslopes: np.ndarray
    downslopes: np.ndarray
    upslope_max_positions: np.ndarray
    downslope_min_positions: np.ndarray


def sloperange(
    decg: ArrayLike,
    r_wave_times: ArrayLike,
    sampling_frequency: float,
) -> SlopeRangeResult:
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
    SlopeRangeResult
        Named result containing five one-dimensional arrays. ``edr`` and the
        two position arrays are aligned with ``r_wave_times``; ``upslopes`` and
        ``downslopes`` are aligned with ``decg``. The result can also be
        unpacked in that field order.

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
    Extrema positions use the zero-based Python and canonical Biosiglib sample
    grid. Beats contribute to both diagnostic slope vectors only when both
    analysis windows are complete. The earliest sample is selected when an
    extreme value is tied.

    Examples
    --------
    >>> import numpy as np
    >>> from biosigpy.ecg import sloperange
    >>> decg = np.zeros(40)
    >>> decg[[9, 10]] = 3.0
    >>> decg[[13, 14]] = -2.0
    >>> result = sloperange(decg, [0.1], 100.0)
    >>> result.edr
    array([5.])
    >>> result.upslope_max_positions
    array([9.])
    >>> edr, upslopes, downslopes, upmaxpos, downminpos = result
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
    upslopes = np.full(derivative_ecg.shape, np.nan, dtype=np.float64)
    downslopes = np.full(derivative_ecg.shape, np.nan, dtype=np.float64)
    upslope_max_positions = np.full(
        times.shape, np.nan, dtype=np.float64
    )
    downslope_min_positions = np.full(
        times.shape, np.nan, dtype=np.float64
    )

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

        upslope_segment = derivative_ecg[upslope_start:upslope_stop]
        downslope_segment = derivative_ecg[
            downslope_start:downslope_stop
        ]
        upslope_max_position = upslope_start + int(
            np.argmax(upslope_segment)
        )
        downslope_min_position = downslope_start + int(
            np.argmin(downslope_segment)
        )

        upslopes[upslope_start:upslope_stop] = upslope_segment
        downslopes[downslope_start:downslope_stop] = downslope_segment
        upslope_max_positions[beat_index] = upslope_max_position
        downslope_min_positions[beat_index] = downslope_min_position
        edr[beat_index] = (
            derivative_ecg[upslope_max_position]
            - derivative_ecg[downslope_min_position]
        )

    return SlopeRangeResult(
        edr=edr,
        upslopes=upslopes,
        downslopes=downslopes,
        upslope_max_positions=upslope_max_positions,
        downslope_min_positions=downslope_min_positions,
    )


def _round_half_away_from_zero(values: np.ndarray) -> np.ndarray:
    return np.copysign(np.floor(np.abs(values) + 0.5), values)
