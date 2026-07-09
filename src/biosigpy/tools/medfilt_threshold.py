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

    x = as_real_vector(x, name="x")
    if x.size == 0:
        raise ValueError("x must not be empty")
    if x.size < 2:
        raise ValueError("x must contain at least two samples")
    if np.any(np.isinf(x)):
        raise ValueError("x must not contain infinite values")

    window = as_integer_scalar(window, name="window")
    if window < 2:
        raise ValueError("window must be greater than or equal to 2")

    factor = as_positive_real_scalar(factor, name="factor")
    max_threshold = as_positive_real_scalar(max_threshold, name="max_threshold")

    window = min(window, x.size)
    half_window = window // 2
    padded = np.concatenate(
        (
            x[:half_window][::-1],
            x,
            x[-half_window:][::-1],
        )
    )

    filter_length = window - 1
    mf = np.empty(x.size, dtype=np.float64)
    for output_index in range(x.size):
        padded_index = output_index + half_window
        start = padded_index - filter_length // 2
        window_values = padded[start : start + filter_length]
        mf[output_index] = np.median(window_values)

    threshold = factor * mf
    threshold[threshold > max_threshold] = max_threshold
    return threshold
