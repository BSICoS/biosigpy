"""Estimate pulse-anchored lower and upper PPG envelopes."""

import warnings

import numpy as np
from scipy.interpolate import PchipInterpolator

from biosigpy.ppg.pulsedel import pulsedel
from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def pulsenv(ppg, sampling_frequency, t_d=None, *, t_a=None, t_b=None,
            window_a=0.25, window_b=0.15, max_gap=np.inf):
    """Interpolate basal and apex amplitudes on the original PPG grid.

    Parameters
    ----------
    ppg : array_like
        PPG samples. A missing sample cannot serve as an anchor.
    sampling_frequency : float
        Positive sampling frequency in Hz. The first PPG sample is at zero seconds.
    t_d : array_like, optional
        Detection times in seconds from the first PPG sample. Supply detections
        or both precomputed fiducial vectors.
    t_a, t_b : array_like, optional
        Apex and basal times in seconds from the first PPG sample. NaNs denote missing points.
    window_a, window_b : float, optional
        Forward and backward search windows in seconds for pulsedel.
    max_gap : float, optional
        Split each envelope before interpolation when consecutive valid
        fiducial times differ by at least this many seconds. Default infinity.

    Returns
    -------
    lower_envelope, upper_envelope : ndarray
        Arrays on the same sample grid as ppg. Each needs two distinct
        anchors per segment. Values outside supported segments are NaN.
    """
    ppg = as_real_vector(np.atleast_1d(ppg), name="ppg")
    if not ppg.size or np.isinf(ppg).any():
        raise ValueError("ppg must be nonempty and contain no infinities")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    last_time = (ppg.size - 1) / fs
    window_a = as_positive_real_scalar(window_a, name="window_a")
    window_b = as_positive_real_scalar(window_b, name="window_b")
    gap = np.asarray(max_gap)
    if gap.ndim != 0 or gap.dtype.kind not in "iuf" or not gap > 0:
        raise ValueError("max_gap must be a positive duration or infinity")
    if t_d is not None and t_a is None and t_b is None:
        detections = _validate_times(t_d, last_time, "t_d")
        if not np.isfinite(detections).any():
            warnings.warn("No pulse detections provided.", stacklevel=2)
            return np.full(ppg.size, np.nan), np.full(ppg.size, np.nan)
        apex, basal = pulsedel(ppg, fs, detections, window_a=window_a,
                               window_b=window_b, points="AB")
    elif t_d is None and t_a is not None and t_b is not None:
        apex = _validate_times(t_a, last_time, "t_a")
        basal = _validate_times(t_b, last_time, "t_b")
    else:
        raise ValueError("supply t_d or both t_a and t_b")

    clock_tolerance = 4 * np.spacing(last_time)
    lower = _interpolate_envelope(ppg, basal, fs, float(gap), clock_tolerance, "lower")
    upper = _interpolate_envelope(ppg, apex, fs, float(gap), clock_tolerance, "upper")
    return lower, upper


def _validate_times(values, last_time, name):
    times = as_real_vector(np.atleast_1d(values), name=name)
    if np.any((times < 0) | (times > last_time)):
        raise ValueError(f"{name} must lie within the signal")
    return times


def _interpolate_envelope(ppg, times, fs, max_gap, clock_tolerance, label):
    """Interpolate each supported group independently using its sampled anchors."""
    envelope = np.full(ppg.size, np.nan)
    times = np.sort(times[np.isfinite(times)])
    anchors, first = np.unique(np.floor(times * fs + 0.5).astype(int), return_index=True)
    valid = np.isfinite(ppg[anchors])
    anchors, times = anchors[valid], times[first[valid]]
    boundaries = np.flatnonzero(np.diff(times) + clock_tolerance >= max_gap) + 1
    for group in np.split(anchors, boundaries):
        if group.size >= 2:
            samples = np.arange(group[0], group[-1] + 1)
            envelope[samples] = PchipInterpolator(group, ppg[group])(samples)
    if not np.isfinite(envelope).any():
        warnings.warn(f"Not enough connected {label} anchors to interpolate.", stacklevel=3)
    return envelope
