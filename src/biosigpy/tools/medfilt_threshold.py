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
    """Compute a capped median-filtered adaptive threshold.

    Parameters
    ----------
    x : array_like
        One-dimensional real numeric signal with at least two samples.
        Infinite values are invalid.
    window : int
        Median window length in samples. Values below 2 are invalid. Windows
        longer than the signal are clipped to the signal length.
    factor : float
        Positive multiplier applied to the median-filtered signal.
    max_threshold : float
        Positive upper bound applied to the threshold.

    Returns
    -------
    numpy.ndarray
        One-dimensional threshold array with the same length as ``x``.

    Raises
    ------
    TypeError
        If numeric inputs are not real numeric data.
    ValueError
        If ``x`` is empty, has fewer than two samples, is not a vector,
        contains infinite values, or if scalar parameters are outside their
        accepted ranges.

    Notes
    -----
    This function implements the Biosiglib ``tools.medfilt_threshold``
    specification.

    Examples
    --------
    >>> import numpy as np
    >>> from biosigpy.tools import medfilt_threshold
    >>> x = np.array([1.0, 1.2, 2.4, 1.1, 1.0])
    >>> medfilt_threshold(x, window=3, factor=1.5, max_threshold=2.0).shape
    (5,)
    """

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
