"""Demonstrate missing-event reconstruction with fillgaps."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.hrv import fillgaps, removefp


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
# Exact output of MATLAB's rng(40), randi(5), and randperm(100, 13).
MATLAB_REMOVED_INDICES_ONE_BASED = np.asarray(
    [6, 29, 30, 44, 51, 55, 59, 63, 64, 73, 79, 90, 93]
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--debug",
        action="store_true",
        help="pause after every interactive reconstruction attempt",
    )
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
    reference = _load_r_wave_times()[:100]
    removed_indices = MATLAB_REMOVED_INDICES_ONE_BASED - 1
    with_gaps = np.delete(reference, removed_indices)
    cleaned = removefp(with_gaps)
    corrected = fillgaps(cleaned, debug=args.debug)

    print("HRV missing-event reconstruction:")
    print("=================================")
    print(f"Reference events:       {reference.size}")
    print(f"Input events:           {with_gaps.size}")
    print(f"After removefp:         {cleaned.size}")
    print(f"Reconstructed events:   {corrected.tn.size}")
    print(f"Unresolved intervals:   {np.count_nonzero(np.isnan(corrected.dtn))}")

    figure, axis = plt.subplots(figsize=(12, 5))
    axis.plot(np.diff(reference), "o-", color="black", label="Reference")
    axis.plot(corrected.dtn, "x-", color="tab:green", label="Filled")
    axis.set_xlabel("Interval index")
    axis.set_ylabel("Interval (s)")
    axis.set_title("Reference and reconstructed RR intervals")
    axis.grid(True)
    axis.legend()
    figure.tight_layout()
    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
