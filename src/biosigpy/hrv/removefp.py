"""False-positive event removal for ordered event-time series."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools import medfilt_threshold
from biosigpy.tools._validation import as_real_vector


__all__ = ["removefp"]


def removefp(tk: ArrayLike) -> np.ndarray:
    """Remove detections preceded by an abnormally short interval.

    Parameters
    ----------
    tk : array_like
        Non-empty, finite, strictly increasing event timestamps in seconds.
        The time origin is unrestricted.

    Returns
    -------
    numpy.ndarray
        One-dimensional event timestamps after simultaneous one-pass removal.

    Raises
    ------
    TypeError
        If ``tk`` is non-numeric or complex.
    ValueError
        If ``tk`` is empty, non-finite, not a vector, or not strictly
        increasing.

    Notes
    -----
    The adaptive baseline uses :func:`biosigpy.tools.medfilt_threshold` with
    the fixed Biosiglib settings. All flags are computed from the original
    interval series before any event is removed. This function does not sort
    its input.

    Examples
    --------
    >>> from biosigpy.hrv import removefp
    >>> removefp([0, 1, 2, 2.2, 3, 4, 5])
    array([0., 1., 2., 3., 4., 5.])
    """

    events = as_real_vector(tk, name="tk")
    if events.size == 0:
        raise ValueError("tk must not be empty")
    if np.any(~np.isfinite(events)):
        raise ValueError("tk must contain only finite values")
    if np.any(np.diff(events) <= 0):
        raise ValueError("tk must be strictly increasing")
    if events.size < 3:
        return events.copy()

    intervals = np.diff(events)
    baseline = medfilt_threshold(
        intervals, window=30, factor=1.0, max_threshold=1.5
    )
    false_positive_intervals = intervals < 0.7 * baseline
    keep = np.concatenate(([True], ~false_positive_intervals))
    return events[keep]
