"""Zero-phase filtering with NaN-aware gap handling."""

import numpy as np
from numpy.typing import ArrayLike
from scipy import signal

from biosigpy.tools._nan_processing import parse_nan_filtering, process_nan_signal


def nan_filtfilt(
    b: ArrayLike, a: ArrayLike, x: ArrayLike, max_gap: int = 0
) -> np.ndarray:
    """Zero-phase filter a signal while preserving long NaN gaps."""

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
