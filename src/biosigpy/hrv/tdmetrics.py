"""Time-domain variability metrics from canonical event times."""

import numpy as np
from numpy.typing import ArrayLike


def tdmetrics(tk: ArrayLike) -> dict[str, float]:
    """Compute time-domain beat or pulse variability metrics.

    Parameters
    ----------
    tk
        Strictly increasing beat or pulse event times in seconds. One-dimensional
        arrays and row or column vectors are accepted.

    Returns
    -------
    dict
        ``mhr``, ``sdnn``, ``sdsd``, ``rmssd``, and ``pnn50``.
    """

    event_times = _real_vector(tk, name="tk")
    if event_times.size == 0:
        raise ValueError("tk must not be empty")
    if not np.all(np.isfinite(event_times)):
        raise ValueError("tk must contain only finite values")
    if np.any(event_times < 0):
        raise ValueError("tk values must be non-negative")

    dtk = np.diff(event_times)
    if np.any(dtk <= 0):
        raise ValueError("tk must be strictly increasing")

    successive_interval_differences = np.diff(dtk)
    return {
        "mhr": float(60.0 / np.mean(dtk)),
        "sdnn": float(1000.0 * np.std(dtk, ddof=1)),
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

    if array.ndim == 1:
        pass
    elif array.ndim == 2 and 1 in array.shape:
        array = array.reshape(-1)
    else:
        raise ValueError(f"{name} must be a vector")

    if not np.issubdtype(array.dtype, np.number) or np.issubdtype(
        array.dtype, np.complexfloating
    ):
        raise TypeError(f"{name} must be real numeric data")

    return np.asarray(array, dtype=np.float64)
