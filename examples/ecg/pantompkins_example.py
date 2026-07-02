"""R-wave detection from the shared Medicom MTD ECG fixture.

The example opens a Matplotlib figure when an interactive backend is available.
If no plot window appears, the active backend is likely non-interactive; run
with ``--save-figure PATH`` to write the figure to disk instead.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.ecg import pantompkins


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
if str(EXAMPLES_ROOT) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_ROOT))

from _support.figures import add_save_figure_argument, show_or_save_figure


FIXTURE_PATH = (
    REPOSITORY_ROOT
    / "examples"
    / "fixtures"
    / "ecg"
    / "medicom_mtd_ecg_respiration.csv"
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_save_figure_argument(parser)
    return parser.parse_args()


def _load_column(path: Path, column_name: str) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None or column_name not in reader.fieldnames:
            raise RuntimeError(f"Column {column_name!r} not found in {path}")
        return np.asarray([float(row[column_name]) for row in reader])


def main() -> None:
    args = _parse_args()

    # Load the same ECG channel used by the Biosigmat Pan-Tompkins example.
    sampling_frequency = 256.0
    ecg = _load_column(FIXTURE_PATH, "ecg")
    time = np.arange(ecg.size, dtype=np.float64) / sampling_frequency

    outputs = pantompkins(ecg, sampling_frequency)
    r_wave_times = outputs["r_wave_times"]
    r_wave_indices = np.rint(r_wave_times * sampling_frequency).astype(int)
    valid_r_wave_indices = (0 <= r_wave_indices) & (r_wave_indices < ecg.size)
    r_wave_plot_times = r_wave_times[valid_r_wave_indices]
    r_wave_indices = r_wave_indices[valid_r_wave_indices]

    print("Pan-Tompkins R-wave Detection:")
    print("==============================")
    print(f"Samples:                   {ecg.size}")
    print(f"Duration:                  {time[-1]:.2f} s")
    print(f"Sampling frequency:        {sampling_frequency:.0f} Hz")
    print(f"Detected R-waves:          {r_wave_times.size}")
    print("Detected R-wave times in seconds:")
    print(np.array2string(r_wave_times, precision=8, separator=", "))
    print("Processing outputs:")
    for output_name in ("ecg_filtered", "decg", "decg_envelope"):
        output_value = outputs[output_name]
        print(f"{output_name}: shape={output_value.shape}")

    figure, axes = plt.subplots(4, 1, sharex=True, figsize=(12, 8))

    axes[0].plot(time, ecg, color="tab:blue", linewidth=1)
    axes[0].plot(
        r_wave_plot_times,
        ecg[r_wave_indices],
        "ro",
        markerfacecolor="red",
        markersize=4,
    )
    axes[0].set_ylabel("ECG")
    axes[0].set_title("Original ECG Signal with Detected R-waves")
    axes[0].grid(True)

    axes[1].plot(time, outputs["ecg_filtered"], color="tab:green", linewidth=1)
    axes[1].plot(
        r_wave_plot_times,
        outputs["ecg_filtered"][r_wave_indices],
        "ro",
        markerfacecolor="red",
        markersize=4,
    )
    axes[1].set_ylabel("Filtered ECG")
    axes[1].set_title("Bandpass Filtered ECG Signal (5-12 Hz)")
    axes[1].grid(True)

    axes[2].plot(time, outputs["decg"], color="tab:purple", linewidth=1)
    axes[2].set_ylabel("Squared Derivative")
    axes[2].set_title("Squared Derivative of Filtered ECG")
    axes[2].grid(True)

    axes[3].plot(time, outputs["decg_envelope"], color="tab:red", linewidth=1)
    axes[3].plot(
        r_wave_plot_times,
        outputs["decg_envelope"][r_wave_indices],
        "ko",
        markerfacecolor="black",
        markersize=4,
    )
    axes[3].set_ylabel("Envelope")
    axes[3].set_xlabel("Time (s)")
    axes[3].set_title("Integrated Envelope with Detected Peaks")
    axes[3].grid(True)

    figure.suptitle("Pan-Tompkins Algorithm: Signal Processing Steps")
    figure.tight_layout()

    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
