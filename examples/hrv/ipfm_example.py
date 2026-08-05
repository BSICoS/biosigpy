"""Demonstrate IPFM heart-rate and TVIPFM modulation reconstruction."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.hrv import ipfm


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
    event_times = _load_r_wave_times()[:100]
    sampling_frequency = 4.0
    result = ipfm(event_times, sampling_frequency, return_m=True)
    sample_times = (
        event_times[0]
        + np.arange(result.ihr.size, dtype=np.float64) / sampling_frequency
    )

    print("IPFM heart-timing reconstruction:")
    print("=================================")
    print(f"Input events:           {event_times.size}")
    print(f"Output samples:         {result.ihr.size}")
    print(f"Sampling frequency:     {sampling_frequency:.1f} Hz")
    print(f"Mean heart rate:        {np.mean(result.ihr):.4f} Hz")

    figure, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].plot(sample_times, result.ihr, color="tab:blue")
    axes[0].set_ylabel("Heart rate (Hz)")
    axes[0].set_title("IPFM-based instantaneous heart rate")

    axes[1].plot(sample_times, result.m, color="tab:red")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Modulating signal")
    axes[1].set_title("TVIPFM normalized modulating signal")

    for axis in axes:
        axis.grid(True)
    figure.tight_layout()
    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
