"""Median-filtered adaptive threshold."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import (
    as_integer_scalar,
    as_positive_real_scalar,
    as_real_vector,
)


def medfilt_threshold(
    x: ArrayLike, window: int, factor: float, max_threshold: float
) -> np.ndarray:
    """Compute a capped median-filtered adaptive threshold."""

    signal = as_real_vector(x, name="x")
    if signal.size == 0:
        raise ValueError("x must not be empty")
    if signal.size < 2:
        raise ValueError("x must contain at least two samples")
    if np.any(np.isinf(signal)):
        raise ValueError("x must not contain infinite values")

    requested_window = as_integer_scalar(window, name="window")
    if requested_window < 2:
        raise ValueError("window must be greater than or equal to 2")

    scale = as_positive_real_scalar(factor, name="factor")
    cap = as_positive_real_scalar(max_threshold, name="max_threshold")

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
