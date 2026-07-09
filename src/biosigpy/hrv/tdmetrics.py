"""Time-domain variability metrics from interval series."""

import numpy as np
from numpy.typing import ArrayLike


def tdmetrics(dtk: ArrayLike) -> dict[str, float]:
    """Compute time-domain beat or pulse variability metrics.

    Parameters
    ----------
    dtk : array_like
        Beat or pulse interval series in seconds. The input must be a
        one-dimensional real numeric vector. ``NaN`` values are treated as
        missing interval markers and are omitted before computing the metrics.
        Infinite, zero, and negative values are invalid.

    Returns
    -------
    dict[str, float]
        Dictionary containing ``mhr`` in beats/min, ``sdnn`` in ms, ``sdsd``
        in ms, ``rmssd`` in ms, and ``pnn50`` in percent.

    Raises
    ------
    TypeError
        If ``dtk`` is not real numeric data.
    ValueError
        If ``dtk`` is empty, is not one-dimensional, contains infinite,
        zero, or negative values, or has fewer than two finite intervals.

    Notes
    -----
    This function implements the Biosiglib ``hrv.tdmetrics`` specification.
    It is modality-generic and can be used with beat or pulse interval series
    that satisfy the input contract.

    Examples
    --------
    >>> import numpy as np
    >>> from biosigpy.hrv import tdmetrics
    >>> metrics = tdmetrics(np.array([0.82, 0.80, 0.84, np.nan, 0.81]))
    >>> sorted(metrics)
    ['mhr', 'pnn50', 'rmssd', 'sdnn', 'sdsd']
    """

    intervals = _real_vector(dtk, name="dtk")
    if intervals.size == 0:
        raise ValueError("dtk must not be empty")
    if np.any(np.isinf(intervals)):
        raise ValueError("dtk must not contain infinite values")
    if np.any(intervals[~np.isnan(intervals)] <= 0):
        raise ValueError("dtk values must be positive or NaN")

    valid_intervals = intervals[~np.isnan(intervals)]
    if valid_intervals.size < 2:
        raise ValueError("dtk must contain at least two valid intervals")

    successive_interval_differences = np.diff(valid_intervals)
    return {
        "mhr": float(60.0 / np.mean(valid_intervals)),
        "sdnn": float(1000.0 * np.std(valid_intervals, ddof=1)),
        "sdsd": float(
            1000.0 * np.std(successive_interval_differences, ddof=1)
        ),
        "rmssd": float(
            1000.0 * np.sqrt(np.mean(successive_interval_differences**2))
        ),
        "pnn50": float(
            100.0
            * np.count_nonzero(np.abs(successive_interval_differences) > 0.05)
            / successive_interval_differences.size
        ),
    }


def _real_vector(values: ArrayLike, *, name: str) -> np.ndarray:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be numeric")

    try:
        array = np.asarray(values)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be numeric") from error

    if array.ndim != 1:
        raise ValueError(f"{name} must be a vector")

    if not np.issubdtype(array.dtype, np.number) or np.issubdtype(
        array.dtype, np.complexfloating
    ):
        raise TypeError(f"{name} must be real numeric data")

    return np.asarray(array, dtype=np.float64)
