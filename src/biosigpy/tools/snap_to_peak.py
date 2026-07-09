"""Snap approximate detections to local signal maxima."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def snap_to_peak(
    ecg: ArrayLike, detections: ArrayLike, window_size: float = 20
) -> np.ndarray:
    """Refine one-based detection positions to local ECG maxima."""

    ecg = as_real_vector(ecg, name="ecg")
    if ecg.size == 0:
        raise ValueError("ecg must not be empty")
    if ecg.size < 2:
        raise ValueError("ecg must contain at least two samples")
    if np.any(np.isinf(ecg)):
        raise ValueError("ecg must not contain infinite values")

    detections = as_real_vector(detections, name="detections")
    if np.any(np.isinf(detections)):
        raise ValueError("detections must not contain infinite values")

    finite_detections = detections[np.isfinite(detections)]
    if np.any((finite_detections < 1) | (finite_detections > ecg.size)):
        raise ValueError("detections must be between 1 and len(ecg)")

    window_size = as_positive_real_scalar(window_size, name="window_size")
    effective_window_size = _round_positive_half_up(window_size)

    refined = np.empty(detections.size, dtype=np.float64)
    for index, detection in enumerate(detections):
        if np.isnan(detection):
            refined[index] = np.nan
            continue

        sample_index = _round_positive_half_up(detection) - 1
        if np.isnan(ecg[sample_index]):
            refined[index] = np.nan
            continue

        segment_start, segment_end = _finite_segment_bounds(ecg, sample_index)
        window_start = max(segment_start, sample_index - effective_window_size)
        window_end = min(segment_end, sample_index + effective_window_size)
        window = ecg[window_start : window_end + 1]
        refined[index] = window_start + int(np.argmax(window)) + 1

    return refined


def _round_positive_half_up(value: float) -> int:
    return int(np.floor(float(value) + 0.5))


def _finite_segment_bounds(ecg: np.ndarray, sample_index: int) -> tuple[int, int]:
    segment_start = sample_index
    while segment_start > 0 and np.isfinite(ecg[segment_start - 1]):
        segment_start -= 1

    segment_end = sample_index
    while segment_end + 1 < ecg.size and np.isfinite(ecg[segment_end + 1]):
        segment_end += 1

    return segment_start, segment_end
