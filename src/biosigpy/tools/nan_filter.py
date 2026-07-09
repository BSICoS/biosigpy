"""Causal filtering with NaN-aware gap handling."""

import numpy as np
from numpy.typing import ArrayLike
from scipy import signal

from biosigpy.tools._nan_processing import parse_nan_filtering, process_nan_signal


def nan_filter(
    b: ArrayLike, a: ArrayLike, x: ArrayLike, max_gap: int = 0
) -> np.ndarray:
    """Filter a one-dimensional signal while preserving long NaN gaps."""

    b, a, x, max_gap = parse_nan_filtering(b, a, x, max_gap)
    minimum_segment_length = max(b.size, a.size)
    return process_nan_signal(
        b, a, x, max_gap, signal.lfilter, minimum_segment_length
    )
