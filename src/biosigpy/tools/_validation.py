"""Internal validation helpers for biosigpy.tools."""

import numpy as np
from numpy.typing import ArrayLike


def as_real_vector(values: ArrayLike, *, name: str) -> np.ndarray:
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


def as_integer_scalar(value: int, *, name: str) -> int:
    if isinstance(value, (str, bytes, bool)):
        raise TypeError(f"{name} must be an integer")

    array = np.asarray(value)
    if array.ndim != 0:
        raise ValueError(f"{name} must be scalar")
    if np.issubdtype(array.dtype, np.bool_):
        raise TypeError(f"{name} must be an integer")
    if not np.issubdtype(array.dtype, np.integer):
        raise TypeError(f"{name} must be an integer")

    return int(array)


def as_positive_real_scalar(value: float, *, name: str) -> float:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be numeric")

    array = np.asarray(value)
    if array.ndim != 0:
        raise ValueError(f"{name} must be scalar")
    if not np.issubdtype(array.dtype, np.number) or np.issubdtype(
        array.dtype, np.complexfloating
    ):
        raise TypeError(f"{name} must be real numeric data")

    scalar = float(array)
    if not np.isfinite(scalar) or scalar <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return scalar