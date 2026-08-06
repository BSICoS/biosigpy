"""Frequency-domain heart-rate variability metrics."""

from __future__ import annotations

import warnings
from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from biosigpy.tools._validation import as_real_vector


__all__ = [
    "FdMetricsResult",
    "FdMetricsWarning",
    "SeparatedFdMetricsResult",
    "fdmetrics",
]


class FdMetricsResult(NamedTuple):
    """Conventional single-spectrum frequency-domain metrics."""

    hf: float
    lf: float
    lfn: float
    lfhf: float


class SeparatedFdMetricsResult(NamedTuple):
    """Frequency-domain metrics from OSP-separated spectra."""

    urlf: float
    re: float
    r: float


class FdMetricsWarning(UserWarning):
    """Structured Biosiglib diagnostic emitted by :func:`fdmetrics`.

    Attributes
    ----------
    warning_id : str
        Canonical Biosiglib warning identifier.
    affected_ids : tuple of str
        Complete set of affected canonical input or output identifiers.
    """

    def __init__(self, warning_id: str, affected_ids: tuple[str, ...]) -> None:
        self.warning_id = warning_id
        self.affected_ids = affected_ids
        affected = ", ".join(affected_ids)
        super().__init__(f"{warning_id}; affected_ids: {affected}")


def fdmetrics(
    pxx: ArrayLike | None = None,
    f: ArrayLike | None = None,
    limit_hf: bool = True,
    *,
    related_pxx: ArrayLike | None = None,
    unrelated_pxx: ArrayLike | None = None,
) -> FdMetricsResult | SeparatedFdMetricsResult:
    """Calculate frequency-domain HRV metrics on a supplied frequency grid.

    Use ``pxx`` and ``f`` for conventional LF/HF metrics. For spectra produced
    after respiration-related OSP decomposition, omit ``pxx`` and provide
    ``related_pxx`` and ``unrelated_pxx`` together with ``f``.

    Parameters
    ----------
    pxx : array_like, optional
        Nonnegative single-spectrum power spectral density. NaN values produce
        an all-NaN single-spectrum result.
    f : array_like
        Finite, nonnegative, strictly increasing frequency samples in hertz.
    limit_hf : bool, default=True
        Limit the high-frequency band to the first sample at or above 0.4 Hz.
        This option applies only to single-spectrum mode.
    related_pxx : array_like, optional
        Nonnegative respiration-related OSP spectrum.
    unrelated_pxx : array_like, optional
        Nonnegative respiration-unrelated OSP spectrum.

    Returns
    -------
    FdMetricsResult or SeparatedFdMetricsResult
        Named, unpackable metrics for the selected call mode.

    Warns
    -----
    FdMetricsWarning
        Emits the canonical ``excessive_vlf_power`` or
        ``zero_required_power`` diagnostic. The warning object exposes its
        canonical ``warning_id`` and complete ``affected_ids`` tuple.

    Raises
    ------
    TypeError
        If numeric inputs or ``limit_hf`` have invalid types.
    ValueError
        If call modes are mixed, vectors have invalid values or shapes, or
        spectrum lengths do not match the frequency grid.

    Examples
    --------
    >>> result = fdmetrics([1, 1, 1], [0.04, 0.15, 0.4])
    >>> result.lf, result.hf
    (0.10999999999999999, 0.25)
    >>> separated = fdmetrics(
    ...     f=[0.04, 0.15, 0.4],
    ...     related_pxx=[0.01, 0.01, 0.01],
    ...     unrelated_pxx=[0.001, 0.001, 0.001],
    ... )
    >>> round(separated.r, 6)
    0.02965
    """

    if not isinstance(limit_hf, (bool, np.bool_)):
        raise TypeError("limit_hf must be a boolean")

    single_mode = (
        pxx is not None and related_pxx is None and unrelated_pxx is None
    )
    separated_mode = (
        pxx is None and related_pxx is not None and unrelated_pxx is not None
    )
    if not single_mode and not separated_mode:
        raise ValueError(
            "provide either pxx or both related_pxx and unrelated_pxx"
        )
    if separated_mode and not bool(limit_hf):
        raise ValueError("limit_hf is not available in separated mode")

    frequencies = as_real_vector(f, name="f")
    _validate_frequencies(frequencies)

    if single_mode:
        spectrum = _validate_spectrum(pxx, name="pxx")
        _require_matching_length(spectrum, frequencies, name="pxx")
        if np.any(np.isnan(spectrum)):
            return _nan_single_result()
        return _single_metrics(spectrum, frequencies, bool(limit_hf))

    related = _validate_spectrum(related_pxx, name="related_pxx")
    unrelated = _validate_spectrum(unrelated_pxx, name="unrelated_pxx")
    _require_matching_length(related, frequencies, name="related_pxx")
    _require_matching_length(unrelated, frequencies, name="unrelated_pxx")
    if np.any(np.isnan(related)) or np.any(np.isnan(unrelated)):
        return _nan_separated_result()
    return _separated_metrics(related, unrelated, frequencies)


def _validate_spectrum(values: ArrayLike | None, *, name: str) -> NDArray[np.float64]:
    spectrum = as_real_vector(values, name=name)
    if spectrum.size == 0:
        raise ValueError(f"{name} must not be empty")
    finite_values = spectrum[~np.isnan(spectrum)]
    if np.any(~np.isfinite(finite_values)) or np.any(finite_values < 0):
        raise ValueError(f"{name} must contain nonnegative finite values or NaN")
    return spectrum


def _validate_frequencies(frequencies: NDArray[np.float64]) -> None:
    if frequencies.size == 0:
        raise ValueError("f must not be empty")
    if np.any(~np.isfinite(frequencies)) or np.any(frequencies < 0):
        raise ValueError("f must contain finite, nonnegative frequencies")
    if np.any(np.diff(frequencies) <= 0):
        raise ValueError("f must be strictly increasing")


def _require_matching_length(
    spectrum: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    *,
    name: str,
) -> None:
    if spectrum.size != frequencies.size:
        raise ValueError(f"{name} and f must have the same length")


def _single_metrics(
    spectrum: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    limit_hf: bool,
) -> FdMetricsResult:
    if _has_excessive_vlf_power(spectrum, frequencies):
        _emit_warning("excessive_vlf_power", ("pxx",))

    band_indices = _required_band_indices(frequencies)
    if band_indices is None:
        return _nan_single_result()
    lf_start, lf_end = band_indices

    if limit_hf and frequencies[-1] >= 0.4:
        hf_end = int(np.flatnonzero(frequencies >= 0.4)[0])
    else:
        hf_end = frequencies.size - 1

    lf = _trapezoid(spectrum[lf_start : lf_end + 1], frequencies[lf_start : lf_end + 1])
    hf = _trapezoid(spectrum[lf_end : hf_end + 1], frequencies[lf_end : hf_end + 1])

    affected_ids: list[str] = []
    if hf == 0:
        affected_ids.append("hf")
    if lf == 0:
        affected_ids.append("lf")
    if affected_ids:
        _emit_warning("zero_required_power", tuple(affected_ids))
        return _nan_single_result()

    return FdMetricsResult(
        hf=hf,
        lf=lf,
        lfn=lf / (lf + hf),
        lfhf=lf / hf,
    )


def _separated_metrics(
    related: NDArray[np.float64],
    unrelated: NDArray[np.float64],
    frequencies: NDArray[np.float64],
) -> SeparatedFdMetricsResult:
    affected_ids = tuple(
        name
        for name, spectrum in (
            ("related_pxx", related),
            ("unrelated_pxx", unrelated),
        )
        if _has_excessive_vlf_power(spectrum, frequencies)
    )
    if affected_ids:
        _emit_warning("excessive_vlf_power", affected_ids)

    band_indices = _required_band_indices(frequencies)
    if band_indices is None:
        return _nan_separated_result()
    lf_start, lf_end = band_indices

    raw_re = _trapezoid(related, frequencies)
    raw_urlf = _trapezoid(
        unrelated[lf_start : lf_end + 1],
        frequencies[lf_start : lf_end + 1],
    )
    if raw_urlf == 0:
        _emit_warning("zero_required_power", ("urlf",))
        return _nan_separated_result()

    re = np.nan if raw_re > 0.05 else raw_re
    urlf = np.nan if raw_urlf > 0.003 else raw_urlf
    ratio = urlf / (re + urlf)
    return SeparatedFdMetricsResult(urlf=urlf, re=re, r=ratio)


def _required_band_indices(
    frequencies: NDArray[np.float64],
) -> tuple[int, int] | None:
    lf_start_matches = np.flatnonzero(frequencies >= 0.04)
    lf_end_matches = np.flatnonzero(frequencies >= 0.15)
    if lf_start_matches.size == 0 or lf_end_matches.size == 0:
        return None
    return int(lf_start_matches[0]), int(lf_end_matches[0])


def _has_excessive_vlf_power(
    spectrum: NDArray[np.float64], frequencies: NDArray[np.float64]
) -> bool:
    if not np.any(frequencies < 0.04):
        return False
    boundary_matches = np.flatnonzero(frequencies >= 0.04)
    if boundary_matches.size == 0:
        return False

    boundary = int(boundary_matches[0])
    vlf_power = _trapezoid(spectrum[: boundary + 1], frequencies[: boundary + 1])
    rest_power = _trapezoid(spectrum[boundary:], frequencies[boundary:])
    if rest_power == 0:
        return vlf_power > 0
    return vlf_power / rest_power > 0.05


def _trapezoid(values: NDArray[np.float64], frequencies: NDArray[np.float64]) -> float:
    if values.size < 2:
        return 0.0
    intervals = np.diff(frequencies)
    return float(np.sum(intervals * (values[:-1] + values[1:]) / 2.0))


def _emit_warning(warning_id: str, affected_ids: tuple[str, ...]) -> None:
    warnings.warn(FdMetricsWarning(warning_id, affected_ids), stacklevel=3)


def _nan_single_result() -> FdMetricsResult:
    return FdMetricsResult(hf=np.nan, lf=np.nan, lfn=np.nan, lfhf=np.nan)


def _nan_separated_result() -> SeparatedFdMetricsResult:
    return SeparatedFdMetricsResult(urlf=np.nan, re=np.nan, r=np.nan)
