"""Integral pulse frequency modulation heart-timing reconstruction."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import make_interp_spline
from scipy.signal import butter, filtfilt

from biosigpy.tools._validation import (
    as_integer_scalar,
    as_positive_real_scalar,
    as_real_vector,
)


__all__ = ["IpfmResult", "ipfm"]


class IpfmResult(NamedTuple):
    """Named, unpackable result returned when the modulating signal is requested.

    Attributes
    ----------
    ihr : numpy.ndarray
        Uniformly sampled instantaneous heart rate in hertz.
    m : numpy.ndarray
        Dimensionless TVIPFM modulating signal.
    """

    ihr: np.ndarray
    m: np.ndarray


def ipfm(
    tn: ArrayLike,
    fs: float,
    spline_order: int = 14,
    *,
    return_m: bool = False,
) -> np.ndarray | IpfmResult:
    """Reconstruct instantaneous heart rate from event timestamps.

    Parameters
    ----------
    tn : array_like
        Finite, strictly increasing event timestamps in seconds. At least two
        timestamps are required.
    fs : float
        Positive output sampling frequency in hertz. When ``return_m=True``,
        it must additionally be greater than 0.06 Hz.
    spline_order : int, default=14
        Order of the canonical B-spline interpolation. It must be between 2
        and the number of event timestamps plus 20, inclusive.
    return_m : bool, default=False
        If true, also compute the TVIPFM modulating signal and return an
        :class:`IpfmResult`. Otherwise, return only instantaneous heart rate.

    Returns
    -------
    numpy.ndarray or IpfmResult
        Instantaneous heart rate alone, or the named ``ihr`` and ``m`` arrays
        when the modulating signal is requested.

    Raises
    ------
    TypeError
        If an input has an invalid numeric or boolean type.
    ValueError
        If timestamps, sampling frequency, spline order, output-grid length,
        or a computed numerical result violates the Biosiglib contract.

    Notes
    -----
    The cumulative beat-count spline uses the canonical ``aptknt`` knot
    sequence. The optional fourth-order 0.03 Hz Butterworth trend is applied
    forward and backward with the fixed MATLAB-compatible odd padding of 12
    samples.

    Examples
    --------
    >>> from biosigpy.hrv import ipfm
    >>> ipfm([0, 1, 2, 3], 4).shape
    (13,)
    >>> result = ipfm([0, 1, 2, 3], 4, return_m=True)
    >>> result.ihr.shape == result.m.shape
    True
    """

    events = as_real_vector(tn, name="tn")
    if events.size < 2:
        raise ValueError("tn must contain at least two event timestamps")
    if np.any(~np.isfinite(events)):
        raise ValueError("tn must contain only finite values")
    if np.any(np.diff(events) <= 0):
        raise ValueError("tn must be strictly increasing")

    sampling_frequency = as_positive_real_scalar(fs, name="fs")
    order = as_integer_scalar(spline_order, name="spline_order")
    if order < 2 or order > events.size + 20:
        raise ValueError(
            "spline_order must be between 2 and the number of extended sites"
        )
    return_m = _as_boolean(return_m, name="return_m")
    if return_m and sampling_frequency <= 0.06:
        raise ValueError("fs must be greater than 0.06 Hz when return_m is true")

    extended_events = _extend_events(events)
    knots = _aptknt(extended_events, order)
    cumulative_beats = np.arange(
        1, extended_events.size + 1, dtype=np.float64
    )
    heart_timing_spline = make_interp_spline(
        extended_events,
        cumulative_beats,
        k=order - 1,
        t=knots,
        check_finite=False,
    )
    heart_timing_spline.extrapolate = False

    candidate_count = int(
        np.ceil((events[-1] - events[0]) * sampling_frequency)
    ) + 1
    sample_indices = np.arange(candidate_count, dtype=np.float64)
    sample_times = events[0] + sample_indices / sampling_frequency
    sample_times = sample_times[sample_times <= events[-1]]
    ihr = np.asarray(
        heart_timing_spline.derivative()(sample_times), dtype=np.float64
    )
    if np.any(~np.isfinite(ihr)) or np.any(ihr <= 0):
        raise ValueError(
            "sampled instantaneous heart rate must be finite and positive"
        )

    if not return_m:
        return ihr
    if ihr.size <= 12:
        raise ValueError(
            "at least 13 sampled values are required for the modulating signal"
        )

    numerator, denominator = butter(
        4, 0.06 / sampling_frequency, btype="lowpass"
    )
    mean_ihr = filtfilt(
        numerator,
        denominator,
        ihr,
        method="pad",
        padtype="odd",
        padlen=12,
    )
    if np.any(~np.isfinite(mean_ihr)) or np.any(mean_ihr <= 0):
        raise ValueError(
            "mean instantaneous heart rate must be finite and positive"
        )
    m = (ihr - mean_ihr) / mean_ihr
    if np.any(~np.isfinite(m)):
        raise ValueError("modulating signal must be finite")
    return IpfmResult(ihr=ihr, m=np.asarray(m, dtype=np.float64))


def _extend_events(events: NDArray[np.float64]) -> NDArray[np.float64]:
    intervals = np.diff(events)
    boundary_interval_count = min(8, intervals.size)
    first_interval = float(np.median(intervals[:boundary_interval_count]))
    last_interval = float(np.median(intervals[-boundary_interval_count:]))
    prepend = events[0] - np.arange(10, 0, -1, dtype=np.float64) * first_interval
    append = events[-1] + np.arange(1, 11, dtype=np.float64) * last_interval
    return np.concatenate((prepend, events, append))


def _aptknt(events: NDArray[np.float64], order: int) -> NDArray[np.float64]:
    interior = np.asarray(
        [
            np.mean(events[index + 1 : index + order])
            for index in range(events.size - order)
        ],
        dtype=np.float64,
    )
    return np.concatenate(
        (
            np.repeat(events[0], order),
            interior,
            np.repeat(events[-1], order),
        )
    )


def _as_boolean(value: bool, *, name: str) -> bool:
    if not isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a boolean")
    return bool(value)
