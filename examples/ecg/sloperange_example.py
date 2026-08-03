"""Slope-range ECG-derived respiration from the shared Medicom MTD fixture."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.ecg import sloperange
from biosigpy.tools import lpd_filter, nan_filter


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
    signals = _load_columns(SIGNAL_FIXTURE_PATH, ("ecg", "respiration"))
    timing = _load_columns(TIMING_FIXTURE_PATH, ("r_wave_times",))
    ecg = signals["ecg"]
    respiration = signals["respiration"]
    r_wave_times = timing["r_wave_times"]

    derivative_filter, _ = lpd_filter(
        sampling_frequency, 50.0, order=4
    )
    derivative_ecg = nan_filter(
        derivative_filter, [1.0], ecg, max_gap=0
    )
    slope_range = sloperange(
        derivative_ecg, r_wave_times, sampling_frequency
    )
    time = np.arange(ecg.size, dtype=np.float64) / sampling_frequency
    r_wave_samples = np.floor(
        r_wave_times * sampling_frequency + 0.5
    ).astype(int)

    print("Slope-range ECG-derived respiration:")
    print("=====================================")
    print(f"ECG samples:                {ecg.size}")
    print(f"R waves:                    {r_wave_times.size}")
    print(
        "Finite EDR estimates:       "
        f"{np.count_nonzero(np.isfinite(slope_range.edr))}"
    )
    print(
        "Boundary NaN estimates:     "
        f"{np.count_nonzero(np.isnan(slope_range.edr))}"
    )

    figure, axes = plt.subplots(4, 1, sharex=True, figsize=(12, 9))
    axes[0].plot(time, ecg, color="tab:blue", linewidth=1)
    axes[0].plot(
        r_wave_times,
        ecg[r_wave_samples],
        "ro",
        markerfacecolor="red",
        markersize=3,
    )
    axes[0].set_ylabel("ECG")
    axes[0].set_title("ECG with reference R waves")

    axes[1].plot(time, derivative_ecg, color="tab:purple", linewidth=1)
    axes[1].plot(
        time,
        slope_range.upslopes,
        color="black",
        linewidth=2,
        label="Upslope windows",
    )
    axes[1].plot(
        time,
        slope_range.downslopes,
        color="tab:pink",
        linewidth=2,
        linestyle="--",
        label="Downslope windows",
    )
    finite_upmax = np.isfinite(slope_range.upslope_max_positions)
    finite_downmin = np.isfinite(slope_range.downslope_min_positions)
    upmax_positions = slope_range.upslope_max_positions[finite_upmax].astype(
        int
    )
    downmin_positions = slope_range.downslope_min_positions[
        finite_downmin
    ].astype(int)
    axes[1].plot(
        time[upmax_positions],
        derivative_ecg[upmax_positions],
        "v",
        color="tab:blue",
        label="Upslope maxima",
    )
    axes[1].plot(
        time[downmin_positions],
        derivative_ecg[downmin_positions],
        "^",
        color="tab:red",
        label="Downslope minima",
    )
    axes[1].set_ylabel("Derivative ECG")
    axes[1].set_title("Slope windows and selected extrema")
    axes[1].legend(loc="upper right")

    axes[2].plot(
        r_wave_times, slope_range.edr, color="tab:green", linewidth=1
    )
    axes[2].set_ylabel("EDR")
    axes[2].set_title("Slope-range ECG-derived respiration")

    axes[3].plot(time, respiration, color="tab:orange", linewidth=1)
    axes[3].set_ylabel("Respiration")
    axes[3].set_xlabel("Time (s)")
    axes[3].set_title("Device respiration reference")

    for axis in axes:
        axis.grid(True)
    figure.tight_layout()
    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
