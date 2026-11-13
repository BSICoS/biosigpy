"""Estimate one-sided PSD from complete windows without crossing NaN gaps."""

import warnings
import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools import interpgap, slicesignal, trimnans
from biosigpy.tools._validation import as_integer_scalar, as_positive_real_scalar, as_real_vector


class NanPwelchWarning(UserWarning):
    """Diagnostic for a channel with no usable analysis windows."""

    def __init__(self, warning_id: str, channel: int):
        self.warning_id = warning_id
        self.affected_ids = ("signal",)
        super().__init__(f"{warning_id}: channel {channel} has no complete valid window")


def nanpwelch(
    signal: ArrayLike, window: ArrayLike | int, overlap: int, nfft: int,
    sampling_frequency: float, max_gap: int | None = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | list[np.ndarray]]:
    """Average complete-window PSDs within independent valid stretches.

    Parameters
    ----------
    signal : array_like
        Nonempty real vector or matrix with samples in rows and channels in
        columns. NaNs mark gaps; infinities are invalid.
    window : int or array_like
        Positive integer length for rectangular weights, or an explicit
        nonzero-energy real coefficient vector.
    overlap : int
        Overlap in samples, less than the window length.
    nfft : int
        Positive DFT length. If shorter than a window, fold weighted samples
        modulo nfft; if longer, zero-pad. No detrending is applied.
    sampling_frequency : float
        Finite positive sample rate in Hz.
    max_gap : int or None, optional
        Maximum internal NaN run to linearly interpolate. Zero, None or an
        empty value disables interpolation.

    Returns
    -------
    psd : ndarray
        One-sided average PSD, a vector for one signal or a frequency-by-channel
        matrix. Keep NaN columns for unusable channels if another is usable.
    frequencies : ndarray
        Frequency bins in Hz, or empty if no channel has a usable window.
    window_psd : ndarray or list of ndarray
        Window-by-frequency spectra for one signal, or one such array per
        channel. Each window has equal weight in the average.
    """
    # Parse and validate inputs
    values = np.asarray(signal)
    if values.dtype.kind not in "iuf":
        raise TypeError("signal must be real numeric data")
    if values.ndim not in (1, 2) or values.size == 0:
        raise ValueError("signal must be a nonempty vector or matrix")
    if np.any(np.isinf(values)):
        raise ValueError("signal must not contain infinity")
    vector_input = values.ndim == 1 or 1 in values.shape
    values = np.asarray(values, dtype=float)
    if vector_input:
        values = values.reshape(-1, 1)
    raw_window = np.asarray(window)
    if raw_window.size == 1:
        if raw_window.dtype.kind not in "iuf" or not np.isfinite(raw_window).all():
            raise ValueError("window length must be finite real numeric data")
        length = float(raw_window.reshape(-1)[0])
        if length < 1 or length != np.floor(length):
            raise ValueError("window length must be a positive integer")
        window_length = int(length)
        weights = np.ones(window_length)
    else:
        weights = as_real_vector(window, name="window")
        if weights.size == 0 or not np.all(np.isfinite(weights)) or not np.any(weights):
            raise ValueError("window must have finite coefficients and positive energy")
        window_length = weights.size
    overlap = as_integer_scalar(overlap, name="overlap")
    nfft = as_integer_scalar(nfft, name="nfft")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    if not 0 <= overlap < window_length or nfft < 1:
        raise ValueError("require 0 <= overlap < window length and nfft > 0")
    if max_gap is None or np.asarray(max_gap).size == 0:
        max_gap = 0
    else:
        max_gap = as_integer_scalar(max_gap, name="max_gap")
        if max_gap < 0:
            raise ValueError("max_gap must be nonnegative")
    if values.shape[0] < window_length:
        raise ValueError("window length must not exceed the original signal length")

    # Process each signal and accumulate results
    spectra = []
    for channel in range(values.shape[1]):
        current_signal = values[:, channel]
        trimmed_signal = trimnans(current_signal)
        if trimmed_signal.size == 0:
            spectra.append(np.empty((0, 0)))
            continue
        if trimmed_signal.size < window_length:
            warnings.warn(NanPwelchWarning("signal_too_short", channel), stacklevel=2)
            spectra.append(np.empty((0, 0)))
            continue

        # Interpolate small gaps and find continuous stretches of valid data
        interpolated_signal = interpgap(trimmed_signal, max_gap) if max_gap else trimmed_signal
        changes = np.diff(np.r_[False, ~np.isnan(interpolated_signal), False].astype(np.int8))
        starts = np.flatnonzero(changes == 1)
        ends = np.flatnonzero(changes == -1)

        # Compute periodograms for all complete windows in signal order
        periodograms = []
        for first, last in zip(starts, ends):
            stretch = interpolated_signal[first:last]
            if stretch.size >= window_length:
                windows = slicesignal(stretch, int(window_length), overlap)
                periodograms.append(_window_periodograms(windows, weights, nfft, fs))
        if not periodograms:
            warnings.warn(NanPwelchWarning("no_valid_segments", channel), stacklevel=2)
            spectra.append(np.empty((0, 0)))
        else:
            spectra.append(np.concatenate(periodograms, axis=0))

    # Average windows independently while keeping all channel positions
    if not any(item.size for item in spectra):
        return np.empty(0), np.empty(0), spectra[0] if vector_input else spectra
    frequencies = np.arange(nfft // 2 + 1) * fs / nfft
    psd = np.full((frequencies.size, values.shape[1]), np.nan)
    for channel, periodograms in enumerate(spectra):
        if periodograms.size:
            psd[:, channel] = np.mean(periodograms, axis=0)
    return (psd[:, 0] if vector_input else psd), frequencies, (spectra[0] if vector_input else spectra)


def _window_periodograms(windows, weights, nfft, fs):
    """Match the MATLAB periodogram PSD, including short-transform folding."""
    weighted = windows * weights
    if nfft < weights.size:
        folded = np.zeros((windows.shape[0], nfft))
        for index in range(weights.size):
            folded[:, index % nfft] += weighted[:, index]
        weighted = folded
    transformed = np.fft.rfft(weighted, n=nfft, axis=1)
    psd = np.abs(transformed)**2 / (fs * np.sum(weights**2))
    psd[:, 1:-1 if nfft % 2 == 0 else None] *= 2
    return psd
