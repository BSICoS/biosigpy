"""Interpolate short internal NaN gaps using original neighboring samples."""

import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import CubicSpline, PchipInterpolator

from biosigpy.tools._validation import as_integer_scalar, as_real_vector


def interpgap(signal: ArrayLike, max_gap: int, method: str = "linear") -> np.ndarray:
    """Reconstruct short internal gaps without changing sample alignment.

    Parameters
    ----------
    signal : array_like
        Real vector with finite values or NaN on a regular sampling grid.
    max_gap : int
        Maximum number of consecutive missing samples to reconstruct.
    method : {'linear', 'nearest', 'pchip', 'spline'}, optional
        Interpolation method; defaults to linear.

    Returns
    -------
    interpolated_signal : ndarray
        One-dimensional copy with eligible gaps reconstructed.

    Raises
    ------
    TypeError
        If signal is not real numeric data or max_gap is not an integer.
    ValueError
        If signal is not a vector, contains infinity, max_gap is negative,
        or method is unsupported.

    Notes
    -----
    Long gaps separate independent segments. Each segment uses only its
    original finite samples as support. Boundary gaps are never extrapolated.
    Nearest ties choose the right sample. Spline uses not-a-knot conditions,
    a quadratic with three knots, and a line with two knots.
    """
    values = as_real_vector(signal, name="signal")
    max_gap = as_integer_scalar(max_gap, name="max_gap")
    if max_gap < 0:
        raise ValueError("max_gap must be nonnegative")
    if not isinstance(method, str) or method not in ("linear", "nearest", "pchip", "spline"):
        raise ValueError("method must be linear, nearest, pchip, or spline")
    if np.any(np.isinf(values)):
        raise ValueError("signal must contain finite values or NaN")
    result = values.copy()
    missing = np.isnan(values)
    changes = np.diff(np.r_[False, missing, False].astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    preserved = (starts == 0) | (ends == values.size) | (ends-starts > max_gap)
    segment_starts = np.r_[0, ends[preserved]]
    segment_ends = np.r_[starts[preserved], values.size]
    for first, last in zip(segment_starts, segment_ends):
        segment = values[first:last]
        targets = np.flatnonzero(np.isnan(segment))
        knots = np.flatnonzero(~np.isnan(segment))
        if targets.size == 0 or knots.size < 2:
            continue
        support = segment[knots]
        if method == "nearest":
            right = np.searchsorted(knots, targets)
            left = right - 1
            choose_right = knots[right]-targets <= targets-knots[left]
            filled = support[np.where(choose_right, right, left)]
        elif method == "linear" or knots.size == 2:
            filled = np.interp(targets, knots, support)
        elif method == "pchip":
            filled = PchipInterpolator(knots, support, extrapolate=False)(targets)
        else:
            filled = CubicSpline(knots, support, bc_type="not-a-knot", extrapolate=False)(targets)
        result[first + targets] = filled
    return result
