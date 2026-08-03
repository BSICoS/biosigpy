"""Iterative reconstruction of missing events in ordered time series."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import PchipInterpolator

from biosigpy.tools import medfilt_threshold
from biosigpy.tools._validation import (
    as_positive_real_scalar,
    as_real_vector,
)


__all__ = ["FillGapsResult", "fillgaps"]


class FillGapsResult(NamedTuple):
    """Named, unpackable result of :func:`fillgaps`.

    Attributes
    ----------
    tn : numpy.ndarray
        Corrected event timestamps in seconds. Every original timestamp is
        preserved exactly and reconstructed timestamps may be inserted.
    dtn : numpy.ndarray
        Successive corrected intervals in seconds. An interval spanning an
        unresolved gap is represented by ``NaN``.
    """

    tn: np.ndarray
    dtn: np.ndarray


def fillgaps(
    tk: ArrayLike,
    gap_detection_factor: float = 1.5,
    correction_upper_factor: float = 1.15,
    correction_lower_factor: float = 0.75,
    minimum_interval: float = 0.5,
    max_gap_duration: float = 10.0,
) -> FillGapsResult:
    """Reconstruct missing event timestamps inside locally detected gaps.

    Parameters
    ----------
    tk : array_like
        Non-empty, finite, strictly increasing event timestamps in seconds.
        The series must already have undergone any desired false-positive
        removal.
    gap_detection_factor : float, default=1.5
        Positive factor applied to the local adaptive baseline for detection.
    correction_upper_factor : float, default=1.15
        Upper factor used to accept reconstructed intervals.
    correction_lower_factor : float, default=0.75
        Lower factor used to detect over-insertion.
    minimum_interval : float, default=0.5
        Non-negative absolute interval boundary in seconds.
    max_gap_duration : float, default=10.0
        Positive maximum gap duration attempted for reconstruction, in
        seconds.

    Returns
    -------
    FillGapsResult
        Named result containing corrected timestamps ``tn`` and their aligned
        successive intervals ``dtn``. The result can also be unpacked in that
        order.

    Raises
    ------
    TypeError
        If a numeric input is non-numeric or complex.
    ValueError
        If timestamps, parameter values, or factor relationships violate the
        Biosiglib contract.

    Notes
    -----
    This function never calls :func:`biosigpy.hrv.removefp` implicitly. The
    recommended explicit preprocessing sequence is ``removefp(tk)`` followed
    by ``fillgaps(cleaned_tk)``. Reconstruction has no plotting or other GUI
    side effects.

    Examples
    --------
    >>> from biosigpy.hrv import fillgaps
    >>> result = fillgaps([0, 1, 2, 4, 5, 6])
    >>> result.tn
    array([0., 1., 2., 3., 4., 5., 6.])
    """

    events = _validate_events(tk)
    gap_detection_factor = as_positive_real_scalar(
        gap_detection_factor, name="gap_detection_factor"
    )
    correction_upper_factor = as_positive_real_scalar(
        correction_upper_factor, name="correction_upper_factor"
    )
    correction_lower_factor = as_positive_real_scalar(
        correction_lower_factor, name="correction_lower_factor"
    )
    minimum_interval = _as_nonnegative_real_scalar(
        minimum_interval, name="minimum_interval"
    )
    max_gap_duration = as_positive_real_scalar(
        max_gap_duration, name="max_gap_duration"
    )
    if not (
        correction_lower_factor
        < correction_upper_factor
        <= gap_detection_factor
    ):
        raise ValueError(
            "factors must satisfy 0 < correction_lower_factor < "
            "correction_upper_factor <= gap_detection_factor"
        )

    if events.size < 3:
        return FillGapsResult(events.copy(), np.diff(events))

    corrected = events.copy()
    unresolved_pairs: set[tuple[float, float]] = set()
    next_insertion_count: dict[tuple[float, float], int] = {}
    previous_reconstructions: dict[
        tuple[float, float], NDArray[np.float64]
    ] = {}

    while True:
        intervals = np.diff(corrected)
        baseline = medfilt_threshold(
            intervals, window=30, factor=1.0, max_threshold=1.5
        )
        pairs = [
            (float(corrected[index]), float(corrected[index + 1]))
            for index in range(intervals.size)
        ]
        pair_to_index = {pair: index for index, pair in enumerate(pairs)}
        detected_indices = {
            index
            for index, interval in enumerate(intervals)
            if interval > gap_detection_factor * baseline[index]
            and interval > minimum_interval
        }
        blocked_indices = {
            pair_to_index[pair]
            for pair in unresolved_pairs
            if pair in pair_to_index
        }
        candidate_indices = sorted(detected_indices - blocked_indices)
        if not candidate_indices:
            break

        all_unresolved_indices = detected_indices | blocked_indices
        accepted: dict[tuple[float, float], NDArray[np.float64]] = {}
        newly_unresolved: set[tuple[float, float]] = set()

        for gap_index in candidate_indices:
            pair = pairs[gap_index]
            gap_duration = intervals[gap_index]
            if gap_duration > max_gap_duration:
                newly_unresolved.add(pair)
                continue

            support = _interpolation_support(
                intervals, gap_index, all_unresolved_indices
            )
            if support is None:
                newly_unresolved.add(pair)
                continue

            insertion_count = next_insertion_count.get(pair, 1)
            reconstruction = _reconstruct_intervals(
                support, gap_duration, insertion_count
            )
            lower_boundary = max(
                correction_lower_factor * baseline[gap_index],
                minimum_interval,
            )
            upper_boundary = correction_upper_factor * baseline[gap_index]

            if np.all(reconstruction < lower_boundary):
                previous = previous_reconstructions.get(pair)
                if previous is None:
                    newly_unresolved.add(pair)
                else:
                    accepted[pair] = previous
                continue

            if np.all(reconstruction < upper_boundary):
                accepted[pair] = reconstruction
                continue

            previous_reconstructions[pair] = reconstruction
            next_insertion_count[pair] = insertion_count + 1

        unresolved_pairs.update(newly_unresolved)
        for pair in newly_unresolved | set(accepted):
            next_insertion_count.pop(pair, None)
            previous_reconstructions.pop(pair, None)

        if accepted:
            corrected = _apply_reconstructions(corrected, accepted)

        if not accepted and not newly_unresolved:
            # Every attempted gap advanced to its next insertion count. The
            # loop recomputes the baseline before the next segment-wide pass.
            continue

        active_pairs = set(pairs[index] for index in candidate_indices)
        for pair in set(next_insertion_count) - active_pairs:
            next_insertion_count.pop(pair, None)
            previous_reconstructions.pop(pair, None)

    corrected_intervals = np.diff(corrected)
    corrected_pairs = [
        (float(corrected[index]), float(corrected[index + 1]))
        for index in range(corrected_intervals.size)
    ]
    for index, pair in enumerate(corrected_pairs):
        if pair in unresolved_pairs:
            corrected_intervals[index] = np.nan
    return FillGapsResult(corrected, corrected_intervals)


def _validate_events(tk: ArrayLike) -> NDArray[np.float64]:
    events = as_real_vector(tk, name="tk")
    if events.size == 0:
        raise ValueError("tk must not be empty")
    if np.any(~np.isfinite(events)):
        raise ValueError("tk must contain only finite values")
    if np.any(np.diff(events) <= 0):
        raise ValueError("tk must be strictly increasing")
    return events


def _as_nonnegative_real_scalar(value: float, *, name: str) -> float:
    if isinstance(value, (str, bytes, bool)):
        raise TypeError(f"{name} must be real numeric data")
    array = np.asarray(value)
    if array.ndim != 0:
        raise ValueError(f"{name} must be scalar")
    if not np.issubdtype(array.dtype, np.number) or np.issubdtype(
        array.dtype, np.complexfloating
    ):
        raise TypeError(f"{name} must be real numeric data")
    scalar = float(array)
    if not np.isfinite(scalar) or scalar < 0:
        raise ValueError(f"{name} must be finite and non-negative")
    return scalar


def _interpolation_support(
    intervals: NDArray[np.float64],
    gap_index: int,
    unresolved_indices: set[int],
) -> NDArray[np.float64] | None:
    excluded = unresolved_indices - {gap_index}
    previous: list[float] = []
    for index in range(gap_index - 1, -1, -1):
        if index not in excluded:
            previous.append(float(intervals[index]))
            if len(previous) == 2:
                break
    following: list[float] = []
    for index in range(gap_index + 1, intervals.size):
        if index not in excluded:
            following.append(float(intervals[index]))
            if len(following) == 2:
                break
    if len(previous) < 2 or len(following) < 2:
        return None
    return np.asarray([*reversed(previous), *following], dtype=np.float64)


def _reconstruct_intervals(
    support: NDArray[np.float64],
    gap_duration: float,
    insertion_count: int,
) -> NDArray[np.float64]:
    coordinates = np.asarray(
        [-1.0, 0.0, insertion_count + 2.0, insertion_count + 3.0]
    )
    targets = np.arange(1, insertion_count + 2, dtype=np.float64)
    raw = PchipInterpolator(coordinates, support, extrapolate=False)(targets)
    reconstructed = raw * gap_duration / np.sum(raw)
    reconstructed[-1] = gap_duration - np.sum(reconstructed[:-1])
    return np.asarray(reconstructed, dtype=np.float64)


def _apply_reconstructions(
    events: NDArray[np.float64],
    reconstructions: dict[tuple[float, float], NDArray[np.float64]],
) -> NDArray[np.float64]:
    output = [float(events[0])]
    for index in range(events.size - 1):
        left = float(events[index])
        right = float(events[index + 1])
        reconstruction = reconstructions.get((left, right))
        if reconstruction is not None:
            output.extend(left + np.cumsum(reconstruction[:-1]))
        output.append(right)
    return np.asarray(output, dtype=np.float64)
