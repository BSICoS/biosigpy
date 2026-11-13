"""Find the highest local peak along a signal axis."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_integer_scalar


def localmax(
    signal: ArrayLike, *, axis: int | None = None
) -> tuple[float, float] | tuple[np.ndarray, np.ndarray]:
    """Return the highest local peak value and position per signal.

    Parameters
    ----------
    signal : array_like
        Nonempty real vector or matrix. NaNs are ignored during detection.
    axis : int, optional
        Zero-based search axis. By default use the first non-singleton
        axis, or axis 0 when all dimensions are singleton. For matrices,
        axis 0 searches columns and axis 1 searches rows.

    Returns
    -------
    peak_values, peak_positions : float or ndarray
        Peak values and zero-based positions in the original input along
        the search axis. Both are NaN when no peak exists. A vector yields
        two scalars; a matrix yields vectors along the unsearched axis.

    Notes
    -----
    A peak must be bounded by strictly lower neighbors. A plateau uses its
    middle retained sample, choosing the earlier middle sample for an even
    length. Ties between highest peaks select the first one.
    """
    values = np.asarray(signal)
    if values.ndim not in (1, 2):
        raise ValueError("signal must be a vector or two-dimensional matrix")
    if values.dtype.kind not in "iuf":
        raise TypeError("signal must be real numeric data")
    if values.size == 0:
        raise ValueError("signal must not be empty")
    if axis is None:
        axis = next((i for i, size in enumerate(values.shape) if size > 1), 0)
    else:
        axis = as_integer_scalar(axis, name="axis")
        if not 0 <= axis < values.ndim:
            raise ValueError("axis is outside the signal dimensions")

    if values.ndim == 1:
        return _highest_peak(values)
    peaks = [_highest_peak(vector) for vector in np.moveaxis(values, axis, -1)]
    heights, positions = zip(*peaks)
    return np.asarray(heights), np.asarray(positions)


def _highest_peak(signal: np.ndarray) -> tuple[float, float]:
    original_positions = np.flatnonzero(~np.isnan(signal))
    samples = signal[original_positions]
    if samples.size < 3:
        return np.nan, np.nan

    starts = np.r_[0, np.flatnonzero(samples[1:] != samples[:-1]) + 1]
    heights = samples[starts]
    peaks = np.flatnonzero((heights[1:-1] > heights[:-2]) &
                          (heights[1:-1] > heights[2:])) + 1
    if peaks.size == 0:
        return np.nan, np.nan

    winner = peaks[np.argmax(heights[peaks])]
    center = (starts[winner] + starts[winner + 1] - 1) // 2
    return float(heights[winner]), float(original_positions[center])
