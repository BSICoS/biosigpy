"""ECG baseline removal from fiducial isoelectric samples."""

from __future__ import annotations

import warnings
from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import CubicSpline

from biosigpy.tools._validation import as_integer_scalar, as_real_vector


__all__ = [
    "BaselineRemoveResult",
    "BaselineRemoveWarning",
    "baselineremove",
]


class BaselineRemoveResult(NamedTuple):
    """Named, unpackable result of :func:`baselineremove`."""

    ecg_detrended: np.ndarray
    baseline: np.ndarray


class BaselineRemoveWarning(UserWarning):
    """Structured Biosiglib diagnostic emitted by :func:`baselineremove`."""

    def __init__(self) -> None:
        self.warning_id = "no_valid_fiducial_positions"
        self.affected_ids = ("fiducial_positions",)
        super().__init__(
            "no_valid_fiducial_positions; "
            "affected_ids: fiducial_positions"
        )


def baselineremove(
    ecg: ArrayLike,
    fiducial_positions: ArrayLike,
    offset: int,
    window_size: int = 5,
) -> BaselineRemoveResult:
    """Remove an ECG baseline estimated at supplied fiducial positions.

    Parameters
    ----------
    ecg : array_like
        Non-empty one-dimensional finite real ECG signal.
    fiducial_positions : array_like
        Non-empty positive finite positions on the canonical one-based sample
        grid. Fractional, unordered, and repeated positions are accepted.
    offset : int
        Nonnegative sample offset subtracted before position rounding.
    window_size : int, default=5
        Positive nominal averaging-window size. Even values produce a
        symmetric effective span of ``window_size + 1`` samples away from
        signal boundaries.

    Returns
    -------
    BaselineRemoveResult
        Named result containing the detrended ECG and estimated baseline as
        one-dimensional arrays aligned with ``ecg``.

    Warns
    -----
    BaselineRemoveWarning
        If no adjusted fiducial position lies inside the signal. The original
        ECG and an all-zero baseline are returned.

    Raises
    ------
    TypeError
        If a vector is non-numeric or complex, or a scalar is not an integer.
    ValueError
        If a vector is empty, has an invalid shape or nonfinite values, a
        scalar is outside its allowed range, or exactly one valid fiducial
        position remains.

    Notes
    -----
    Positions are rounded half away from zero after subtracting ``offset``.
    Two valid positions define a line, three define a quadratic, and four or
    more define a cubic not-a-knot spline. The end pieces extrapolate over ECG
    samples outside the first and last valid positions.

    Examples
    --------
    >>> result = baselineremove([2, 4, 8, 16, 32], [1, 5], 0, 2)
    >>> result.baseline
    array([ 3.  ,  8.25, 13.5 , 18.75, 24.  ])
    """

    signal = as_real_vector(ecg, name="ecg")
    if signal.size == 0:
        raise ValueError("ecg must not be empty")
    if not np.all(np.isfinite(signal)):
        raise ValueError("ecg must contain only finite values")

    positions = as_real_vector(
        fiducial_positions, name="fiducial_positions"
    )
    if positions.size == 0:
        raise ValueError("fiducial_positions must not be empty")
    if not np.all(np.isfinite(positions)):
        raise ValueError(
            "fiducial_positions must contain only finite values"
        )
    if np.any(positions <= 0.0):
        raise ValueError("fiducial_positions must be positive")

    offset_value = as_integer_scalar(offset, name="offset")
    if offset_value < 0:
        raise ValueError("offset must be nonnegative")
    window_value = as_integer_scalar(window_size, name="window_size")
    if window_value <= 0:
        raise ValueError("window_size must be positive")

    adjusted_positions = _round_half_away_from_zero(
        positions - offset_value
    )
    valid_positions = np.unique(
        adjusted_positions[
            (adjusted_positions >= 1) &
            (adjusted_positions <= signal.size)
        ]
    ).astype(np.intp)

    if valid_positions.size == 0:
        warnings.warn(BaselineRemoveWarning(), stacklevel=2)
        return BaselineRemoveResult(signal.copy(), np.zeros_like(signal))
    if valid_positions.size == 1:
        raise ValueError("at least two valid fiducial positions are required")

    radius = window_value // 2
    fiducial_levels = np.asarray(
        [
            np.mean(
                signal[
                    max(0, int(position) - 1 - radius):
                    min(signal.size, int(position) + radius)
                ]
            )
            for position in valid_positions
        ],
        dtype=np.float64,
    )
    spline = CubicSpline(
        valid_positions.astype(np.float64),
        fiducial_levels,
        bc_type="not-a-knot",
        extrapolate=True,
    )
    sample_grid = np.arange(1, signal.size + 1, dtype=np.float64)
    baseline = np.asarray(spline(sample_grid), dtype=np.float64)
    return BaselineRemoveResult(signal - baseline, baseline)


def _round_half_away_from_zero(values: np.ndarray) -> np.ndarray:
    return np.copysign(np.floor(np.abs(values) + 0.5), values)
