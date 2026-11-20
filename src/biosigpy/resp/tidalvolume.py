"""Estimate respiratory amplitude and its peak and valley envelopes."""

import numpy as np
from scipy.interpolate import PchipInterpolator

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def tidalvolume(resp, sampling_frequency, mindist=0, *, max_gap=np.inf):
    """Estimate peak-to-following-valley amplitude on the respiratory sample grid.

    Parameters
    ----------
    resp : array_like
        Nonempty, uniformly sampled respiratory signal. NaNs interrupt detection.
    sampling_frequency : float
        Sampling frequency in Hz.
    mindist : float, optional
        Minimum separation in seconds between consecutive raw zero crossings.
        The default of zero accepts every crossing.
    max_gap : float, optional
        Split each interpolation before fitting when successive anchors for
        that curve are at least this many seconds apart. Defaults to infinity.

    Returns
    -------
    tdvol, upper_envelope, lower_envelope : ndarray
        Cycle amplitudes at peak-valley midpoints, peak amplitudes and valley
        amplitudes, each interpolated on the input grid in input amplitude units.
        Samples outside groups of at least two anchors are NaN.
    """
    resp = as_real_vector(np.atleast_1d(resp), name="resp").copy()
    if not resp.size or np.isinf(resp).any():
        raise ValueError("resp must be nonempty and contain no infinities")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    mindist = np.asarray(mindist)
    if mindist.ndim or mindist.dtype.kind not in "iuf" or not np.isfinite(mindist) or mindist < 0:
        raise ValueError("mindist must be a finite nonnegative duration")
    max_gap = np.asarray(max_gap)
    if max_gap.ndim or max_gap.dtype.kind not in "iuf" or not max_gap > 0:
        raise ValueError("max_gap must be a positive duration or infinity")

    # Remove the linear trend using observed samples at their original times.
    valid = np.isfinite(resp)
    observed = np.flatnonzero(valid)
    if observed.size >= 2:
        centered_time = observed - observed.mean()
        centered_resp = resp[valid] - resp[valid].mean()
        slope = (centered_time @ centered_resp) / (centered_time @ centered_time)
        resp[valid] = centered_resp - slope * centered_time

    # Detect complete half-cycles and pair extrema within each observed segment.
    edges = np.diff(np.r_[False, valid, False].astype(int))
    peaks, valleys, volumes = [], [], []
    for start, stop in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        segment = resp[start:stop]
        crossings, directions = _detect_zero_crossings(segment, fs, mindist)
        pending_peak = None
        for first, last, direction in zip(crossings[:-1], crossings[1:], directions):
            samples = segment[first:last + 1]
            if direction > 0:
                index = start + first + np.argmax(samples)
                peaks.append((index, resp[index]))
                pending_peak = index
            else:
                index = start + first + np.argmin(samples)
                valleys.append((index, resp[index]))
                if pending_peak is not None:
                    midpoint = (pending_peak + index) / 2
                    volumes.append((midpoint, resp[pending_peak] - resp[index]))
                    pending_peak = None

    # Fit each curve on the sample grid from its own anchors.
    tdvol = _interpolate_envelope(volumes, resp.size, fs, max_gap)
    upper = _interpolate_envelope(peaks, resp.size, fs, max_gap)
    lower = _interpolate_envelope(valleys, resp.size, fs, max_gap)
    return tdvol, upper, lower


def _detect_zero_crossings(segment, fs, mindist):
    """Locate changes of sign and keep the first of consecutive equal directions."""
    nonzero = np.flatnonzero(segment)
    changes = np.diff(np.sign(segment[nonzero]))
    indices = np.flatnonzero(changes)
    crossings = nonzero[indices + 1]
    directions = changes[indices]
    if crossings.size:
        keep = np.r_[True, np.diff(crossings) / fs >= mindist]
        crossings, directions = crossings[keep], directions[keep]
        alternating = np.r_[True, np.diff(directions) != 0]
        crossings, directions = crossings[alternating], directions[alternating]
    return crossings, directions


def _interpolate_envelope(anchors, size, fs, max_gap):
    """Interpolate each group only between its first and last anchor."""
    envelope = np.full(size, np.nan)
    if len(anchors) < 2:
        return envelope
    anchors = np.asarray(anchors)
    breaks = np.flatnonzero(np.diff(anchors[:, 0]) / fs >= max_gap) + 1
    for group in np.split(anchors, breaks):
        if len(group) >= 2:
            samples = np.arange(int(np.ceil(group[0, 0])), int(np.floor(group[-1, 0])) + 1)
            envelope[samples] = PchipInterpolator(group[:, 0], group[:, 1])(samples)
    return envelope
