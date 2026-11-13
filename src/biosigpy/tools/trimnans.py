"""Remove boundary NaN samples from a signal."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_real_vector


def trimnans(signal: ArrayLike) -> np.ndarray:
    """Remove leading and trailing NaNs, preserving internal gaps.

    Parameters
    ----------
    signal : array_like
        One-dimensional real signal containing finite values or NaN.
        Row and column vectors represent the same sequence.

    Returns
    -------
    trimmed_signal : ndarray
        Independent one-dimensional copy of the inclusive first-to-last
        finite span. Empty and all-NaN inputs return an empty vector.

    Raises
    ------
    TypeError
        If the signal is not real numeric data.
    ValueError
        If the signal is not a vector or contains infinite values.

    Notes
    -----
    No interpolation is performed. Preserve the original time offset
    separately if the trimmed signal must remain aligned with other data.
    """
    values = as_real_vector(signal, name="signal")
    if np.any(np.isinf(values)):
        raise ValueError("signal must contain finite values or NaN")
    finite_positions = np.flatnonzero(~np.isnan(values))
    if finite_positions.size == 0:
        return np.empty(0, dtype=np.float64)
    return values[finite_positions[0] : finite_positions[-1] + 1].copy()
