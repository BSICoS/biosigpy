"""Respiration-related HRV decomposition by orthogonal projection."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


__all__ = ["OspResult", "osp"]


class OspResult(NamedTuple):
    """Named, unpackable output of :func:`osp`."""

    m_resp: np.ndarray | float
    m_unrelated: np.ndarray | float
    delay: int | np.ndarray


def osp(
    m: ArrayLike,
    resp: ArrayLike,
    resp_pxx: ArrayLike,
    f: ArrayLike,
    fs: float,
    min_resp_frequency: float = 0.1,
) -> OspResult:
    """Separate respiration-related and unrelated HRV modulation.

    Parameters
    ----------
    m : array_like
        Uniformly sampled, dimensionless HRV modulating signal.
    resp : array_like
        Respiration samples aligned with ``m`` on the same time grid.
    resp_pxx : array_like
        Finite, nonnegative respiration power spectral density.
    f : array_like
        Finite, strictly increasing frequency samples in hertz.
    fs : float
        Positive sampling frequency in hertz.
    min_resp_frequency : float, default=0.1
        Positive lower bound for the selected respiratory frequency.

    Returns
    -------
    OspResult
        Respiration-related modulation, orthogonal residual, and adaptive
        delayed-respiration model order. The component arrays align with
        ``m[delay - 1:]``.

    Raises
    ------
    TypeError
        If an input has an invalid numeric type.
    ValueError
        If vector shapes, spectrum values, frequencies, sampling parameters,
        or finite signal lengths violate the Biosiglib contract.

    Notes
    -----
    The Gram-matrix pseudoinverse uses the explicit Biosiglib binary64
    threshold, rather than NumPy's default pseudoinverse tolerance.

    Examples
    --------
    >>> result = osp(
    ...     [99, 1, 2, 3, 4, 5],
    ...     [1, 0, -1, 0, 1, 0],
    ...     [0, 0, 1],
    ...     [0, 0.5, 1],
    ...     1,
    ... )
    >>> result.delay
    2
    >>> np.allclose(result.m_resp + result.m_unrelated, [1, 2, 3, 4, 5])
    True
    """

    modulation = as_real_vector(m, name="m")
    respiration = as_real_vector(resp, name="resp")
    spectrum = as_real_vector(resp_pxx, name="resp_pxx")
    frequencies = as_real_vector(f, name="f")
    sampling_frequency = as_positive_real_scalar(fs, name="fs")
    minimum_frequency = as_positive_real_scalar(
        min_resp_frequency, name="min_resp_frequency"
    )

    _validate_spectrum(spectrum, frequencies)

    if modulation.size == 0 or respiration.size == 0:
        empty = np.asarray([], dtype=np.float64)
        return OspResult(empty.copy(), empty.copy(), empty.copy())
    if np.any(np.isnan(modulation)) or np.any(np.isnan(respiration)):
        empty = np.asarray([], dtype=np.float64)
        return OspResult(empty.copy(), empty.copy(), empty.copy())
    if np.any(np.isinf(modulation)) or np.any(np.isinf(respiration)):
        raise ValueError("m and resp must not contain infinite values")
    if modulation.size != respiration.size:
        raise ValueError("m and resp must have the same length")

    low_frequency, high_frequency = _occupied_power_limits(
        spectrum, frequencies
    )
    band_mask = (frequencies >= low_frequency) & (
        frequencies <= high_frequency
    )
    if not np.any(band_mask):
        band_mask = np.ones(frequencies.shape, dtype=bool)

    band_spectrum = spectrum[band_mask]
    band_frequencies = frequencies[band_mask]
    dominant_frequency = _dominant_frequency(
        band_spectrum, band_frequencies
    )
    dominant_frequency = max(dominant_frequency, minimum_frequency)
    delay = max(
        int(np.floor(2.0 * sampling_frequency / dominant_frequency + 0.5)),
        1,
    )

    if modulation.size < delay:
        return OspResult(np.nan, np.nan, delay)

    subspace = np.lib.stride_tricks.sliding_window_view(
        respiration, delay
    )
    gram = subspace.T @ subspace
    gram_inverse = _gram_pseudoinverse(gram)
    projection = subspace @ gram_inverse @ subspace.T
    delayed_modulation = modulation[delay - 1 :]
    m_resp = projection @ delayed_modulation
    m_unrelated = delayed_modulation - m_resp
    return OspResult(
        m_resp=np.asarray(m_resp, dtype=np.float64),
        m_unrelated=np.asarray(m_unrelated, dtype=np.float64),
        delay=delay,
    )


def _validate_spectrum(
    spectrum: NDArray[np.float64], frequencies: NDArray[np.float64]
) -> None:
    if spectrum.size < 2 or frequencies.size < 2:
        raise ValueError("resp_pxx and f must contain at least two samples")
    if np.any(~np.isfinite(spectrum)) or np.any(spectrum < 0):
        raise ValueError("resp_pxx must contain finite, nonnegative values")
    if np.any(~np.isfinite(frequencies)):
        raise ValueError("f must contain finite values")
    if spectrum.size != frequencies.size:
        raise ValueError("resp_pxx and f must have the same length")
    if np.any(np.diff(frequencies) <= 0):
        raise ValueError("f must be strictly increasing")


def _occupied_power_limits(
    spectrum: NDArray[np.float64], frequencies: NDArray[np.float64]
) -> tuple[float, float]:
    average_delta = (frequencies[-1] - frequencies[0]) / (
        frequencies.size - 1
    )
    widths = np.empty(frequencies.shape, dtype=np.float64)
    if frequencies[0] == 0:
        widths[:-1] = np.diff(frequencies)
        widths[-1] = average_delta
    else:
        widths[0] = average_delta
        widths[1:] = np.diff(frequencies)

    cumulative_power = np.concatenate(
        ([0.0], np.cumsum(spectrum * widths, dtype=np.float64))
    )
    boundaries = np.concatenate(
        (
            [frequencies[0]],
            (frequencies[:-1] + frequencies[1:]) / 2.0,
            [frequencies[-1]],
        )
    )
    total_power = cumulative_power[-1]
    if total_power == 0:
        return np.nan, np.nan

    return (
        _interpolate_power_threshold(
            0.05 * total_power, cumulative_power, boundaries
        ),
        _interpolate_power_threshold(
            0.95 * total_power, cumulative_power, boundaries
        ),
    )


def _interpolate_power_threshold(
    threshold: float,
    cumulative_power: NDArray[np.float64],
    boundaries: NDArray[np.float64],
) -> float:
    index = int(np.flatnonzero(threshold <= cumulative_power)[0])
    index = max(index, 1)
    fraction = (threshold - cumulative_power[index - 1]) / (
        cumulative_power[index] - cumulative_power[index - 1]
    )
    return float(
        boundaries[index - 1]
        + fraction * (boundaries[index] - boundaries[index - 1])
    )


def _dominant_frequency(
    spectrum: NDArray[np.float64], frequencies: NDArray[np.float64]
) -> float:
    if spectrum.size < 3:
        return float(frequencies[int(np.argmax(spectrum))])

    peak_indices = _peak_indices(spectrum)
    if not peak_indices:
        return float(frequencies[int(np.argmax(spectrum))])
    if len(peak_indices) <= 3:
        peak_values = spectrum[peak_indices]
        selected = peak_indices[int(np.argmax(peak_values))]
        return float(frequencies[selected])
    return float(frequencies[peak_indices[0]])


def _peak_indices(spectrum: NDArray[np.float64]) -> list[int]:
    peaks: list[int] = []
    index = 1
    while index < spectrum.size - 1:
        if spectrum[index] <= spectrum[index - 1]:
            index += 1
            continue

        plateau_end = index
        while (
            plateau_end + 1 < spectrum.size
            and spectrum[plateau_end + 1] == spectrum[index]
        ):
            plateau_end += 1
        if (
            plateau_end < spectrum.size - 1
            and spectrum[plateau_end] > spectrum[plateau_end + 1]
        ):
            peaks.append(index)
        index = plateau_end + 1
    return peaks


def _gram_pseudoinverse(gram: NDArray[np.float64]) -> NDArray[np.float64]:
    left, singular_values, right_transpose = np.linalg.svd(
        gram, full_matrices=False
    )
    tolerance = max(gram.shape) * np.spacing(singular_values[0])
    reciprocal = np.zeros(singular_values.shape, dtype=np.float64)
    retained = singular_values > tolerance
    reciprocal[retained] = 1.0 / singular_values[retained]
    return (right_transpose.T * reciprocal) @ left.T
