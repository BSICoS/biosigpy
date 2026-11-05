"""Extract fixed-length windows from a sampled signal."""

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_integer_scalar, as_positive_real_scalar, as_real_vector


def slicesignal(
    signal: ArrayLike,
    window_length: int,
    overlap: int,
    sampling_frequency: float | None = None,
    *,
    use_last: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Slice a signal into rows with optional overlapping samples.

    Parameters
    ----------
    signal : array_like
        Nonempty real vector. Existing NaNs and infinities are preserved.
    window_length : int
        Positive window length in samples.
    overlap : int
        Number of shared samples, from zero to window_length minus one.
    sampling_frequency : float, optional
        Finite positive sampling frequency in Hz. If supplied, also return
        one time center per window, in seconds from the first input sample.
    use_last : bool, optional
        Pad an uncovered tail on the right with NaNs. A nonempty signal
        shorter than one window is then returned as one padded window.

    Returns
    -------
    windows : ndarray
        Array of shape (number of windows, window_length).
    time_centers : ndarray, optional
        Returned together with windows only when sampling_frequency is
        supplied. Each center refers to the nominal window, including padding.
    """
    values = as_real_vector(signal, name="signal")
    window_length = as_integer_scalar(window_length, name="window_length")
    overlap = as_integer_scalar(overlap, name="overlap")
    if window_length <= 0 or not 0 <= overlap < window_length:
        raise ValueError("require 0 <= overlap < window_length")
    if not isinstance(use_last, (bool, np.bool_)):
        raise TypeError("use_last must be boolean")
    if sampling_frequency is not None:
        sampling_frequency = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    if values.size == 0:
        raise ValueError("signal must not be empty")
    if values.size < window_length and not use_last:
        raise ValueError("signal is too short for a complete window")

    step = window_length - overlap
    complete_steps, remainder = divmod(max(0, values.size - window_length), step)
    count = 1 + complete_steps + (use_last and remainder != 0)
    padding = (count - 1) * step + window_length - values.size
    if padding > 0:
        values = np.pad(values, (0, padding), constant_values=np.nan)
    windows = sliding_window_view(values, window_length)[::step].copy()
    if sampling_frequency is None:
        return windows
    centers = (np.arange(count) * step + (window_length - 1) / 2) / sampling_frequency
    return windows, centers
