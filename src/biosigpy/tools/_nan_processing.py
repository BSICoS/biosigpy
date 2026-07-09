"""Shared NaN-gap processing for filtering tools."""

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_integer_scalar, as_real_vector


FilterFunc = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]


def parse_nan_filtering(
    b: ArrayLike, a: ArrayLike, x: ArrayLike, max_gap: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    b = as_real_vector(b, name="b")
    a = as_real_vector(a, name="a")
    x = as_real_vector(x, name="x")
    max_gap = as_integer_scalar(max_gap, name="max_gap")

    if b.size == 0:
        raise ValueError("b must contain at least one coefficient")
    if a.size == 0:
        raise ValueError("a must contain at least one coefficient")
    if np.any(~np.isfinite(b)):
        raise ValueError("b must contain finite values")
    if np.any(~np.isfinite(a)):
        raise ValueError("a must contain finite values")
    if np.any(np.isinf(x)):
        raise ValueError("x must not contain infinite values")
    if max_gap < 0:
        raise ValueError("max_gap must be greater than or equal to 0")

    return b, a, x, max_gap


def process_nan_signal(
    b: np.ndarray,
    a: np.ndarray,
    x: np.ndarray,
    max_gap: int,
    filter_func: FilterFunc,
    minimum_segment_length: int,
) -> np.ndarray:
    y = np.full(x.shape, np.nan, dtype=np.float64)
    if x.size == 0:
        return y

    idx_nan = np.isnan(x)
    if np.all(idx_nan):
        return y

    preserved_nan_seqs = _preserved_nan_sequences(idx_nan, max_gap)
    valid_segments = _valid_segments(x.size, preserved_nan_seqs)

    for segment_start, segment_end in valid_segments:
        segment_data = x[segment_start : segment_end + 1]
        segment_filled = _fill_short_gaps(segment_data)

        if np.any(np.isnan(segment_filled)):
            continue
        if segment_filled.size < minimum_segment_length:
            continue

        y[segment_start : segment_end + 1] = filter_func(b, a, segment_filled)

    return y


def _preserved_nan_sequences(idx_nan: np.ndarray, max_gap: int) -> list[tuple[int, int]]:
    preserved_nan_seqs: list[tuple[int, int]] = []
    for segment_start, segment_end in _nan_sequences(idx_nan):
        gap_length = segment_end - segment_start + 1
        is_boundary_gap = segment_start == 0 or segment_end == idx_nan.size - 1
        is_long_internal_gap = not is_boundary_gap and gap_length > max_gap
        if is_boundary_gap or is_long_internal_gap:
            preserved_nan_seqs.append((segment_start, segment_end))
    return preserved_nan_seqs


def _nan_sequences(idx_nan: np.ndarray) -> list[tuple[int, int]]:
    seqs: list[tuple[int, int]] = []
    sample_index = 0
    while sample_index < idx_nan.size:
        if not idx_nan[sample_index]:
            sample_index += 1
            continue

        segment_start = sample_index
        while sample_index + 1 < idx_nan.size and idx_nan[sample_index + 1]:
            sample_index += 1
        seqs.append((segment_start, sample_index))
        sample_index += 1

    return seqs


def _valid_segments(
    signal_length: int, preserved_nan_seqs: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    valid_segments: list[tuple[int, int]] = []
    segment_start = 0
    for gap_start, gap_end in preserved_nan_seqs:
        segment_end = gap_start - 1
        if segment_start <= segment_end:
            valid_segments.append((segment_start, segment_end))
        segment_start = gap_end + 1

    if segment_start < signal_length:
        valid_segments.append((segment_start, signal_length - 1))

    return valid_segments


def _fill_short_gaps(segment_data: np.ndarray) -> np.ndarray:
    if not np.any(np.isnan(segment_data)):
        return segment_data

    sample_indices = np.arange(segment_data.size)
    finite_indices = sample_indices[np.isfinite(segment_data)]
    if (
        finite_indices.size < 2
        or finite_indices[0] != 0
        or finite_indices[-1] != segment_data.size - 1
    ):
        return segment_data.copy()

    return np.interp(
        sample_indices,
        finite_indices,
        segment_data[finite_indices],
    )
