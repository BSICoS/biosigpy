"""Demonstrate fixture-based frequency-domain HRV metrics."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import detrend, welch

from biosigpy.hrv import fdmetrics, ipfm, osp


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = REPOSITORY_ROOT / "examples" / "fixtures" / "ecg"
TIMING_FIXTURE_PATH = FIXTURE_ROOT / "medicom_mtd_r_wave_timing.csv"
RESPIRATION_FIXTURE_PATH = FIXTURE_ROOT / "medicom_mtd_ecg_respiration.csv"


def _load_columns(path: Path, *columns: str) -> tuple[np.ndarray, ...]:
    with path.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None or not set(columns) <= set(reader.fieldnames):
            missing = sorted(set(columns) - set(reader.fieldnames or ()))
            raise RuntimeError(f"Columns {missing} not found in {path}")
        rows = list(reader)
    return tuple(
        np.asarray([float(row[column]) for row in rows], dtype=np.float64)
        for column in columns
    )


def _power_spectrum(
    signal: np.ndarray, sampling_frequency: float
) -> tuple[np.ndarray, np.ndarray]:
    window_length = min(256, signal.size)
    return welch(
        signal,
        fs=sampling_frequency,
        window=np.hamming(window_length),
        nperseg=window_length,
        noverlap=window_length // 2,
        detrend=False,
    )


def main() -> None:
    (event_times,) = _load_columns(TIMING_FIXTURE_PATH, "r_wave_times")
    respiration_times, respiration_signal = _load_columns(
        RESPIRATION_FIXTURE_PATH, "time", "respiration"
    )

    event_times = event_times[:100]
    sampling_frequency = 4.0
    modulation = ipfm(
        event_times, sampling_frequency, return_m=True
    ).m
    sample_times = event_times[0] + (
        np.arange(modulation.size, dtype=np.float64) / sampling_frequency
    )
    aligned_respiration = PchipInterpolator(
        respiration_times,
        detrend(respiration_signal),
        extrapolate=False,
    )(sample_times)

    frequencies, spectrum = _power_spectrum(modulation, sampling_frequency)
    respiration_frequencies, respiration_spectrum = _power_spectrum(
        aligned_respiration, sampling_frequency
    )
    if not np.array_equal(respiration_frequencies, frequencies):
        raise RuntimeError("HRV and respiration spectra must share a frequency grid")

    classic_metrics = fdmetrics(spectrum, frequencies)
    unlimited_metrics = fdmetrics(spectrum, frequencies, limit_hf=False)

    decomposition = osp(
        modulation,
        aligned_respiration,
        respiration_spectrum,
        frequencies,
        sampling_frequency,
    )
    osp_frequencies, related_spectrum = _power_spectrum(
        decomposition.m_resp, sampling_frequency
    )
    unrelated_frequencies, unrelated_spectrum = _power_spectrum(
        decomposition.m_unrelated, sampling_frequency
    )
    if not np.array_equal(unrelated_frequencies, osp_frequencies):
        raise RuntimeError("OSP component spectra must share a frequency grid")
    osp_metrics = fdmetrics(
        f=osp_frequencies,
        related_pxx=related_spectrum,
        unrelated_pxx=unrelated_spectrum,
    )

    print("Classic LF/HF metrics:")
    print(f"  LF   = {classic_metrics.lf:.4f}")
    print(f"  HF   = {classic_metrics.hf:.4f}")
    print(f"  LFn  = {classic_metrics.lfn:.4f}")
    print(f"  LFHF = {classic_metrics.lfhf:.4f}\n")

    print("Unlimited-HF LF/HF metrics:")
    print(f"  LF   = {unlimited_metrics.lf:.4f}")
    print(f"  HF   = {unlimited_metrics.hf:.4f}")
    print(f"  LFn  = {unlimited_metrics.lfn:.4f}")
    print(f"  LFHF = {unlimited_metrics.lfhf:.4f}\n")

    print("OSP-based metrics:")
    print(f"  UrLF = {osp_metrics.urlf:.4f}")
    print(f"  Re   = {osp_metrics.re:.4f}")
    print(f"  R    = {osp_metrics.r:.4f}")


if __name__ == "__main__":
    main()
