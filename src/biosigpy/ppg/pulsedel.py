"""Delineate PPG pulse fiducials relative to detection times."""

import numpy as np
from scipy.signal import find_peaks

from biosigpy.tools import localmax
from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def pulsedel(ppg, sampling_frequency, t_d, dppg=None, *, window_a=0.25, window_b=0.15,
             points="ABMO"):
    """Return apex, basal, mid-amplitude and onset times in seconds.

    Parameters
    ----------
    ppg : array_like
        PPG samples, starting at time zero. NaNs delimit valid spans.
    sampling_frequency : float
        Sampling frequency in Hz.
    t_d : array_like
        Detection times in seconds, in input order, including repeats and NaNs.
    dppg : array_like, optional
        Derivative aligned sample by sample with the PPG. The caller chooses
        the differentiator and compensates any delay.
    window_a, window_b : float, optional
        Apex search duration after the detection and basal search duration
        before it, in seconds.

    points : str, optional
        Distinct letters A, B, M and O, case-insensitive. Return the requested
        points in the given order. One letter returns an array. Default "ABMO".

    Returns
    -------
    t_a, t_b, t_m, t_o : ndarray
        By default, four one-dimensional arrays, one entry per detection. The midpoint
        uses half the amplitude between the sampled basal point and apex.
        Omitting dppg silently returns an all-NaN t_o of the same length.
        An empty t_d gives an empty array for each requested point.
    """
    ppg = as_real_vector(np.atleast_1d(ppg), name="ppg")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    t_d = as_real_vector(np.atleast_1d(t_d), name="t_d")
    window_a = as_positive_real_scalar(window_a, name="window_a")
    window_b = as_positive_real_scalar(window_b, name="window_b")
    if not ppg.size or np.isinf(ppg).any():
        raise ValueError("ppg must be nonempty and contain no infinities")
    if np.any((t_d < 0) | (t_d > (ppg.size - 1) / fs)):
        raise ValueError("detection times must lie within the signal")
    if dppg is not None:
        dppg = as_real_vector(np.atleast_1d(dppg), name="dppg")
        if dppg.size != ppg.size or np.isinf(dppg).any():
            raise ValueError("dppg must match the PPG length and contain no infinities")

    if not isinstance(points, str):
        raise TypeError("points must be a string")
    points = points.upper()
    if not points or set(points) - set("ABMO") or len(set(points)) != len(points):
        raise ValueError("points must contain distinct letters from ABMO")
    outputs = {key: np.full(t_d.size, np.nan) for key in points}
    need_apex = "A" in points or "M" in points or ("O" in points and dppg is not None)
    need_basal = "B" in points or "M" in points
    forward = int(np.floor(min(window_a * fs, ppg.size) + 0.5))
    backward = int(np.floor(min(window_b * fs, ppg.size) + 0.5))
    for k in np.flatnonzero(~np.isnan(t_d)):
        detection = int(np.floor(t_d[k] * fs + 0.5))
        if np.isnan(ppg[detection]):
            continue

        # Anchor both search windows to this detection's finite PPG span.
        basal = None
        if need_basal:
            first = max(0, detection - backward)
            gaps = np.flatnonzero(np.isnan(ppg[first:detection]))
            if gaps.size:
                first += gaps[-1] + 1
            _, offset = localmax(-ppg[first:detection + 1])
            if np.isfinite(offset):
                basal = first + int(offset)
                if "B" in outputs:
                    outputs["B"][k] = _refine_time(ppg, basal, fs)
        if not need_apex:
            continue
        last = min(ppg.size - 1, detection + forward)
        gaps = np.flatnonzero(np.isnan(ppg[detection:last + 1]))
        if gaps.size:
            last = detection + gaps[0] - 1
        _, apex_offset = localmax(ppg[detection:last + 1])
        if np.isnan(apex_offset):
            continue
        apex = detection + int(apex_offset)
        if "A" in outputs:
            outputs["A"][k] = _refine_time(ppg, apex, fs)

        if "M" in outputs and basal is not None:
            target = (ppg[basal] + ppg[apex]) / 2
            distance = np.abs(ppg[basal:apex + 1] - target)
            midpoint = np.argmin(distance)
            outputs["M"][k] = basal / fs + _refine_time(-distance, midpoint, fs)
        if "O" in outputs and dppg is not None:
            outputs["O"][k] = _find_onset(ppg, dppg, apex, fs)

    if len(points) == 1:
        return outputs[points]
    return tuple(outputs.values())


def _refine_time(signal, sample, fs):
    """Refine a sampled extremum using its two neighbors."""
    if sample == 0 or sample == signal.size - 1:
        return sample / fs
    left, center, right = signal[sample - 1:sample + 2]
    curvature = left - 2 * center + right
    if curvature == 0:
        return np.nan
    return (sample + (left - right) / (2 * curvature)) / fs


def _find_onset(ppg, dppg, apex, fs):
    """Find the 5 percent upslope fiducial in the span ending at the apex."""
    first = max(0, apex - int(np.floor(min(0.3 * fs, ppg.size) + 0.5)))
    window = slice(first, apex + 1)
    gaps = np.flatnonzero(np.isnan(ppg[window]) | np.isnan(dppg[window]))
    if gaps.size:
        first += gaps[-1] + 1
    if first >= apex:
        return np.nan

    rising = dppg[first:apex + 1]
    upslope = np.argmax(rising)
    if rising[upslope] <= 0:
        return np.nan
    rising = rising[:upslope + 1]
    threshold = 0.05 * rising[-1]
    if np.any(rising <= threshold):
        onset = np.argmin(np.abs(rising - threshold))
    else:
        minima, _ = find_peaks(-rising)
        onset = minima[-1] if minima.size else np.argmin(rising)
    return (first + onset) / fs
