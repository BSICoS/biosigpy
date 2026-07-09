"""Causal filtering with NaN-aware gap handling."""

import numpy as np
from numpy.typing import ArrayLike
from scipy import signal

from biosigpy.tools._nan_processing import parse_nan_filtering, process_nan_signal


def nan_filter(
    b: ArrayLike, a: ArrayLike, x: ArrayLike, max_gap: int = 0
) -> np.ndarray:
    """Causally filter a signal with NaN-aware gap handling.

    Parameters
    ----------
    b : array_like
        Numerator filter coefficients. Values must be finite.
    a : array_like
        Denominator filter coefficients. Values must be finite.
    x : array_like
        One-dimensional real numeric signal. ``NaN`` values are allowed.
        Infinite values are invalid.
    max_gap : int, optional
        Maximum internal ``NaN`` gap length to interpolate before filtering.

    Returns
    -------
    numpy.ndarray
        One-dimensional filtered signal with the same length as ``x``.

    Raises
    ------
    TypeError
        If ``b``, ``a``, or ``x`` is not real numeric data.
    ValueError
        If filter coefficients are empty or non-finite, if ``x`` contains
        infinite values, or if ``max_gap`` is negative.

    Notes
    -----
    This function implements the Biosiglib ``tools.nan_filter`` specification.
    Filtering is causal and uses :func:`scipy.signal.lfilter` on each valid
    finite segment.

    Boundary ``NaN`` gaps are preserved. Internal short gaps with length less
    than or equal to ``max_gap`` are interpolated before filtering. Internal
    long gaps with length greater than ``max_gap`` are preserved and split the
    signal into separate segments. Segments that are too short for causal
    filtering return ``NaN``.

    Examples
    --------
    >>> import numpy as np
    >>> from biosigpy.tools import nan_filter
    >>> b = np.array([0.5, 0.5])
    >>> a = np.array([1.0])
    >>> x = np.array([1.0, np.nan, 3.0, 4.0])
    >>> nan_filter(b, a, x, max_gap=1).shape
    (4,)
    """

    b, a, x, max_gap = parse_nan_filtering(b, a, x, max_gap)
    minimum_segment_length = max(b.size, a.size)
    return process_nan_signal(
        b, a, x, max_gap, signal.lfilter, minimum_segment_length
    )
