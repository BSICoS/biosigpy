"""Median-filtered adaptive threshold."""

import numpy as np
from numpy.typing import ArrayLike


def medfilt_threshold(
    x: ArrayLike, window: int, factor: float, max_threshold: float
) -> np.ndarray:
    """Compute a capped median-filtered adaptive threshold."""

    signal = _real_vector(x, name="x")
    if signal.size == 0:
        raise ValueError("x must not be empty")
    if signal.size < 2:
        raise ValueError("x must contain at least two samples")
    if np.any(np.isinf(signal)):
        raise ValueError("x must not contain infinite values")

    requested_window = _integer_scalar(window, name="window")
    if requested_window < 2:
        raise ValueError("window must be greater than or equal to 2")

    scale = _positive_scalar(factor, name="factor")
    cap = _positive_scalar(max_threshold, name="max_threshold")

    effective_window = min(requested_window, signal.size)
    half_window = effective_window // 2
    padded = np.concatenate(
        (
            signal[:half_window][::-1],
            signal,
            signal[-half_window:][::-1],
        )
    )

    filter_length = effective_window - 1
    baseline = np.empty(signal.size, dtype=np.float64)
    for output_index in range(signal.size):
        padded_index = output_index + half_window
        start = padded_index - filter_length // 2
        window_values = padded[start : start + filter_length]
        baseline[output_index] = np.median(window_values)

    threshold = scale * baseline
    threshold[threshold > cap] = cap
    return threshold


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


def _integer_scalar(value: int, *, name: str) -> int:
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


def _positive_scalar(value: float, *, name: str) -> float:
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
