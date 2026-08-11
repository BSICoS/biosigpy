"""Baseline-wander removal from the shared Medicom MTD ECG fixture."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.ecg import baselineremove


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
if str(EXAMPLES_ROOT) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_ROOT))

from _support.figures import add_save_figure_argument, show_or_save_figure


SIGNAL_FIXTURE_PATH = (
    REPOSITORY_ROOT
    / "examples"
    / "fixtures"
    / "ecg"
    / "medicom_mtd_ecg_respiration.csv"
)
TIMING_FIXTURE_PATH = (
    REPOSITORY_ROOT
    / "examples"
    / "fixtures"
    / "ecg"
    / "medicom_mtd_r_wave_timing.csv"
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_save_figure_argument(parser)
    return parser.parse_args()


def _load_columns(
    path: Path, column_names: tuple[str, ...]
) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None or not set(column_names).issubset(
            reader.fieldnames
        ):
            raise RuntimeError(f"Columns {column_names!r} not found in {path}")
        rows = list(reader)
    return {
        column_name: np.asarray(
            [float(row[column_name]) for row in rows], dtype=np.float64
        )
        for column_name in column_names
    }


def main() -> None:
    args = _parse_args()
    sampling_frequency = 256.0
    signals = _load_columns(SIGNAL_FIXTURE_PATH, ("time", "ecg"))
    timing = _load_columns(
        TIMING_FIXTURE_PATH, ("r_wave_times", "r_wave_samples")
    )
    time = signals["time"]
    ecg = signals["ecg"]
    r_wave_times = timing["r_wave_times"]
    r_wave_samples = timing["r_wave_samples"].astype(np.intp)

    # Match the Biosigmat example: PR fiducials 150 ms before each R wave.
    offset = int(np.floor(0.15 * sampling_frequency + 0.5))
    result = baselineremove(ecg, r_wave_samples, offset)
    fiducial_samples = r_wave_samples - offset

    # Convert the canonical one-based samples only for NumPy array indexing.
    r_wave_indices = r_wave_samples - 1
    fiducial_indices = fiducial_samples - 1

    print("ECG Baseline-Wander Removal:")
    print("============================")
    print(f"ECG samples:                {ecg.size}")
    print(f"Sampling frequency:         {sampling_frequency:.0f} Hz")
    print(f"R waves:                    {r_wave_samples.size}")
    print(f"PR fiducial offset:         {offset} samples")

    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 7))

    axes[0].plot(time, ecg, color="tab:blue", linewidth=1.5, label="ECG Signal")
    axes[0].plot(
        time,
        result.baseline,
        "r--",
        linewidth=2,
        label="Estimated Baseline",
    )
    axes[0].plot(
        r_wave_times,
        ecg[r_wave_indices],
        "go",
        markerfacecolor="green",
        markersize=4,
        label="R-peaks",
    )
    axes[0].plot(
        time[fiducial_indices],
        ecg[fiducial_indices],
        "ko",
        markerfacecolor="black",
        markersize=4,
        label="Fiducial Points",
    )
    axes[0].set_title("Original ECG Signal with Estimated Baseline")
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Amplitude")
    axes[0].legend(loc="best")
    axes[0].grid(True)

    axes[1].plot(
        time,
        result.ecg_detrended,
        color="tab:blue",
        linewidth=1.5,
        label="Corrected ECG",
    )
    axes[1].plot(
        r_wave_times,
        result.ecg_detrended[r_wave_indices],
        "go",
        markerfacecolor="green",
        markersize=4,
        label="R-peaks",
    )
    axes[1].plot(
        time[fiducial_indices],
        result.ecg_detrended[fiducial_indices],
        "ko",
        markerfacecolor="black",
        markersize=4,
        label="Fiducial Points",
    )
    axes[1].set_title("Baseline-Corrected ECG Signal")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Amplitude")
    axes[1].legend(loc="best")
    axes[1].grid(True)

    figure.tight_layout()
    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
