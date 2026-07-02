"""Time-domain HRV metrics from the shared Medicom MTD timing fixture."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from biosigpy.hrv import tdmetrics


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_PATH = (
    REPOSITORY_ROOT
    / "examples"
    / "fixtures"
    / "ecg"
    / "medicom_mtd_r_wave_timing.csv"
)


def _load_column(path: Path, column_name: str) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None or column_name not in reader.fieldnames:
            raise RuntimeError(f"Column {column_name!r} not found in {path}")
        return np.asarray([float(row[column_name]) for row in reader])


def _medfilt_threshold(
    values: np.ndarray, window: int, factor: float, max_threshold: float
) -> np.ndarray:
    if values.size < window:
        window = values.size

    half_window = window // 2
    padded_values = np.concatenate(
        (
            np.flip(values[:half_window]),
            values,
            np.flip(values[-half_window:]),
        )
    )

    # Match the MATLAB example's medfiltThreshold(dtk, 50, 1.5, 1.5)
    # preprocessing without adding a public Biosigpy helper.
    filter_size = window - 1
    filter_half_width = filter_size // 2
    threshold = np.empty(values.size, dtype=np.float64)
    for index, center in enumerate(range(half_window, half_window + values.size)):
        start = center - filter_half_width
        stop = center + filter_half_width + 1
        threshold[index] = np.median(padded_values[start:stop])
    return np.minimum(factor * threshold, max_threshold)


def main() -> None:
    # Load pre-calculated ECG R-wave timing data from the copied fixture.
    r_wave_times = _load_column(FIXTURE_PATH, "r_wave_times")

    # Compute RR intervals, remove median-threshold outliers, and pass the
    # cleaned interval series to tdmetrics, matching the Biosigmat example.
    intervals = np.diff(r_wave_times)
    threshold = _medfilt_threshold(intervals, 50, 1.5, 1.5)
    outliers = intervals > threshold
    intervals_without_outliers = intervals[~outliers]

    metrics = tdmetrics(intervals_without_outliers)

    print("Time Domain HRV Metrics:")
    print("========================")
    print(f"Mean Heart Rate (mhr):    {metrics['mhr']:.2f} beats/min")
    print(f"SDNN:                     {metrics['sdnn']:.2f} ms")
    print(f"SDSD:                     {metrics['sdsd']:.2f} ms")
    print(f"RMSSD:                    {metrics['rmssd']:.2f} ms")
    print(f"pNN50:                    {metrics['pnn50']:.2f} %")


if __name__ == "__main__":
    main()
