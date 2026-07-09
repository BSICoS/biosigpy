"""Zero-phase filtering with NaN-aware gap handling."""

import numpy as np
from numpy.typing import ArrayLike
from scipy import signal

from biosigpy.tools._nan_processing import parse_nan_filtering, process_nan_signal


def nan_filtfilt(
    b: ArrayLike, a: ArrayLike, x: ArrayLike, max_gap: int = 0
) -> np.ndarray:
    """Zero-phase filter a signal with NaN-aware gap handling.

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
        One-dimensional zero-phase filtered signal with the same length as
        ``x``.

    Raises
    ------
    TypeError
        If ``b``, ``a``, or ``x`` is not real numeric data.
    ValueError
        If filter coefficients are empty or non-finite, if ``x`` contains
        infinite values, or if ``max_gap`` is negative.

    Notes
    -----
    This function implements the Biosiglib ``tools.nan_filtfilt``
    specification. Filtering is zero-phase and uses
    :func:`scipy.signal.filtfilt` on each valid finite segment.

    Boundary ``NaN`` gaps are preserved. Internal short gaps with length less
    than or equal to ``max_gap`` are interpolated before filtering. Internal
    long gaps with length greater than ``max_gap`` are preserved and split the
    signal into separate segments. Segments that are too short for zero-phase
    filtering return ``NaN``.

    Examples
    --------
    >>> import numpy as np
    >>> from biosigpy.tools import nan_filtfilt
    >>> b = np.array([0.25, 0.5, 0.25])
    >>> a = np.array([1.0])
    >>> x = np.linspace(0.0, 1.0, 20)
    >>> x[8] = np.nan
    >>> nan_filtfilt(b, a, x, max_gap=1).shape
    (20,)
    """

    b, a, x, max_gap = parse_nan_filtering(b, a, x, max_gap)
    filter_order = max(b.size - 1, a.size - 1)
    minimum_segment_length = 3 * filter_order + 1
    padlen = 3 * filter_order

    def filter_func(
        b: np.ndarray, a: np.ndarray, segment_filled: np.ndarray
    ) -> np.ndarray:
        return signal.filtfilt(b, a, segment_filled, padlen=padlen)

    return process_nan_signal(
        b, a, x, max_gap, filter_func, minimum_segment_length
    )
