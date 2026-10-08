"""Expand original NaN gaps by a symmetric time radius."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def expandnans(signal: ArrayLike, sampling_frequency: float, expansion_duration: float) -> np.ndarray:
    """Expand original NaN gaps without changing sample alignment.

    Parameters
    ----------
    signal : array_like
        One-dimensional real signal with finite values or NaN.
    sampling_frequency : float
        Finite positive sampling frequency in Hz.
    expansion_duration : float
        Finite nonnegative margin in seconds, applied on each side.

    Returns
    -------
    expanded_signal : ndarray
        One-dimensional copy with expanded NaN spans and unchanged length.

    Raises
    ------
    TypeError
        If any input is not real numeric data.
    ValueError
        If the signal is not a vector, contains infinity, or a scalar
        input is nonfinite or outside its permitted range.

    Notes
    -----
    Exact half-way sample radii round upward. Only original NaNs generate
    expansion; overlapping spans are united and clipped at the boundaries.
    """
    values = as_real_vector(signal, name="signal")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    raw_duration = np.asarray(expansion_duration)
    if not np.issubdtype(raw_duration.dtype, np.number) or np.iscomplexobj(raw_duration):
        raise TypeError("expansion_duration must be real numeric data")
    if raw_duration.ndim != 0:
        raise ValueError("expansion_duration must be scalar")
    duration = float(raw_duration)
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("expansion_duration must be finite and nonnegative")
    if np.any(np.isinf(values)):
        raise ValueError("signal must contain finite values or NaN")
    result = values.copy()
    original_nan = np.isnan(values)
    if not np.any(original_nan):
        return result
    radius = int(np.floor(min(fs * duration, values.size) + 0.5))
    changes = np.diff(np.r_[False, original_nan, False].astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    for start, end in zip(starts, ends):
        result[max(0, start-radius) : min(values.size, end+radius)] = np.nan
    return result
