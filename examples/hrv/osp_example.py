"""Demonstrate respiration-related HRV decomposition with OSP."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import detrend, welch

from biosigpy.hrv import ipfm, osp


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
if str(EXAMPLES_ROOT) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_ROOT))

from _support.figures import add_save_figure_argument, show_or_save_figure


FIXTURE_ROOT = REPOSITORY_ROOT / "examples" / "fixtures" / "ecg"
TIMING_FIXTURE_PATH = FIXTURE_ROOT / "medicom_mtd_r_wave_timing.csv"
RESPIRATION_FIXTURE_PATH = FIXTURE_ROOT / "medicom_mtd_ecg_respiration.csv"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_save_figure_argument(parser)
    return parser.parse_args()


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


def main() -> None:
    args = _parse_args()
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
    window_length = min(256, aligned_respiration.size)
    frequencies, respiration_spectrum = welch(
        aligned_respiration,
        fs=sampling_frequency,
        window=np.hamming(window_length),
        nperseg=window_length,
        noverlap=window_length // 2,
        detrend=False,
    )

    result = osp(
        modulation,
        aligned_respiration,
        respiration_spectrum,
        frequencies,
        sampling_frequency,
    )
    component_times = sample_times[result.delay - 1 :]
    reconstruction_error = np.max(
        np.abs(result.m_resp + result.m_unrelated - modulation[result.delay - 1 :])
    )

    print("Respiration-related HRV decomposition:")
    print("======================================")
    print(f"Adaptive delay:          {result.delay} samples")
    print(f"Component samples:       {result.m_resp.size}")
    print(f"Reconstruction error:    {reconstruction_error:.3e}")

    figure, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True)
    axes[0].plot(
        sample_times,
        aligned_respiration,
        color=(0.24, 0.35, 0.74),
        linewidth=1.4,
    )
    axes[0].set_ylabel(r"$r(n)$")
    axes[0].set_title("OSP Decomposition of HRV Modulating Signal")

    axes[1].plot(
        sample_times,
        modulation,
        color=(0.70, 0.70, 0.72),
        linewidth=2.2,
    )
    axes[1].set_ylabel(r"$m(n)$")

    axes[2].plot(
        component_times,
        result.m_unrelated,
        color=(0.20, 0.20, 0.20),
        linewidth=1.5,
    )
    axes[2].set_ylabel(r"$\hat{m}_{\perp}(n)$")

    axes[3].plot(
        component_times,
        result.m_resp,
        color=(0.20, 0.20, 0.20),
        linewidth=1.5,
        linestyle=":",
    )
    axes[3].set_ylabel(r"$\hat{m}_{r}(n)$")
    axes[3].set_xlabel("Time (seconds)")
    for axis in axes:
        axis.grid(True)
    figure.tight_layout()
    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
