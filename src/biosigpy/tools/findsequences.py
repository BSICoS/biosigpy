"""Find repeated adjacent values with zero-based inclusive boundaries."""

import numpy as np
from numpy.typing import ArrayLike


def findsequences(signal: ArrayLike) -> np.ndarray:
    """Return maximal runs as rows [value, start, end, length].

    Parameters
    ----------
    signal : array_like
        Real numeric vector. NaN and both signs of infinity are accepted.

    Returns
    -------
    sequences : ndarray
        Float array of shape (n, 4). Start and end indices are zero-based
        and inclusive. No runs produces an array of shape (0, 4).

    Raises
    ------
    TypeError
        If the input is not real numeric data.
    ValueError
        If the input is not a vector.

    Notes
    -----
    Runs contain at least two exactly equal adjacent samples. Consecutive
    NaNs compare equal for detection.

    Rows are grouped in this order: nonzero finite values, zeros, NaNs,
    +Inf, and -Inf. Within each group, rows are ordered by start index.

    A result row selects ``signal[int(start):int(end) + 1]``. Equality is
    evaluated before conversion to the floating-point output representation.
    """
    if isinstance(signal, (str, bytes)):
        raise TypeError("signal must be real numeric data")
    values = np.asarray(signal)
    if values.ndim == 2 and 1 in values.shape:
        values = values.reshape(-1)
    if values.ndim != 1:
        raise ValueError("signal must be a vector")
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("signal must be real numeric data")
    if values.size < 2:
        return np.empty((0, 4), dtype=float)
    equal = (values[:-1] == values[1:]) | (np.isnan(values[:-1]) & np.isnan(values[1:]))
    changes = np.diff(np.r_[False, equal, False].astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    repeated = values[starts]
    categories = [np.isfinite(repeated) & (repeated != 0), repeated == 0,
                  np.isnan(repeated), np.isposinf(repeated), np.isneginf(repeated)]
    order = np.concatenate([np.flatnonzero(mask) for mask in categories])
    return np.column_stack((repeated[order], starts[order], ends[order],
                            ends[order] - starts[order] + 1)).astype(float)
