"""Compute Hjorth descriptors from signal and difference energies."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def hjorth(signal: ArrayLike, sampling_frequency: float) -> tuple[float, float, float]:
    """Return mean-square activity, mobility, and frequency spread.

    Parameters
    ----------
    signal : array_like
        Nonempty real vector with finite values or NaN. The baseline is
        not removed. A NaN sample makes all three descriptors undefined.
    sampling_frequency : float
        Finite positive sampling frequency in Hz.

    Returns
    -------
    activity : float
        Mean of the squared input samples.
    mobility : float
        Square root of the first-difference to signal energy ratio,
        scaled by sampling_frequency/(2*pi), or NaN if undefined.
    complexity : float
        Frequency-spread descriptor from the first- and second-difference
        energy ratios, scaled by sampling_frequency/(2*pi), or NaN if
        undefined. It is not a dimensionless mobility ratio.
    """
    values = as_real_vector(signal, name="signal")
    fs = as_positive_real_scalar(sampling_frequency, name="sampling_frequency")
    if values.size == 0 or np.any(np.isinf(values)):
        raise ValueError("signal must be nonempty and contain finite values or NaN")
    # Compute first and second differences
    differences = np.diff(values)
    second_differences = np.diff(differences)

    # Compute signal and difference energies
    e0 = np.sum(values**2)
    e1 = np.sum(differences**2)
    e2 = np.sum(second_differences**2)

    # Compute Hjorth parameters
    activity = float(e0 / values.size)  # Mean-square power
    mobility = np.nan
    complexity = np.nan
    if values.size >= 2 and e0 > 0:
        mobility = float(np.sqrt(e1 / e0) * fs / (2 * np.pi))
    if values.size >= 3 and e0 > 0 and e1 > 0:
        complexity = float(np.sqrt(abs(e2 / e1 - e1 / e0)) * fs / (2 * np.pi))

    return activity, mobility, complexity
