"""Demonstrate false-positive event removal."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.hrv import removefp


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
if str(EXAMPLES_ROOT) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_ROOT))

from _support.figures import add_save_figure_argument, show_or_save_figure


TIMING_FIXTURE_PATH = (
    REPOSITORY_ROOT
    / "examples"
    / "fixtures"
    / "ecg"
    / "medicom_mtd_r_wave_timing.csv"
)
MATLAB_FP_INDICES_ONE_BASED = np.asarray([10, 20, 30])
FP_OFFSETS = np.asarray([0.05, 0.08, 0.06])


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_save_figure_argument(parser)
    return parser.parse_args()


def _load_r_wave_times() -> np.ndarray:
    with TIMING_FIXTURE_PATH.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        if (
            reader.fieldnames is None
            or "r_wave_times" not in reader.fieldnames
        ):
            raise RuntimeError(
                f"Column 'r_wave_times' not found in {TIMING_FIXTURE_PATH}"
            )
        return np.asarray(
            [float(row["r_wave_times"]) for row in reader],
            dtype=np.float64,
        )


def main() -> None:
    args = _parse_args()
    reference = _load_r_wave_times()[:50]
    fp_indices = MATLAB_FP_INDICES_ONE_BASED - 1
    false_positives = reference[fp_indices] + FP_OFFSETS
    observed = np.sort(np.concatenate((reference, false_positives)))

    cleaned = removefp(observed)

    print("HRV false-positive removal:")
    print("===========================")
    print(f"Reference events:       {reference.size}")
    print(f"Observed events:        {observed.size}")
    print(f"After removefp:         {cleaned.size}")

    figure, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=False)
    axes[0].eventplot(
        [reference, observed, cleaned],
        lineoffsets=[3, 2, 1],
        linelengths=0.7,
        colors=["black", "tab:red", "tab:green"],
    )
    axes[0].set_yticks([1, 2, 3])
    axes[0].set_yticklabels(["Cleaned", "Observed", "Reference"])
    axes[0].set_xlabel("Time (s)")
    axes[0].set_title("Reference, observed, and cleaned events")

    axes[1].plot(np.diff(observed), "o-", color="tab:red", label="Observed")
    axes[1].plot(np.diff(cleaned), "o-", color="tab:green", label="Cleaned")
    axes[1].set_xlabel("Interval index")
    axes[1].set_ylabel("Interval (s)")
    axes[1].set_title("Intervals before and after false-positive removal")
    axes[1].legend()
    for axis in axes:
        axis.grid(True)
    figure.tight_layout()
    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
