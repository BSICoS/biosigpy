"""Visualize false-positive removal and missing-event reconstruction."""

from __future__ import annotations

import argparse
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


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_save_figure_argument(parser)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    reference = np.arange(0.0, 24.8, 0.8)
    observed = np.delete(reference, [10, 11, 22])
    observed = np.sort(np.append(observed, reference[5] + 0.12))

    cleaned = removefp(observed)
    corrected = fillgaps(cleaned)

    print("HRV event preprocessing:")
    print("========================")
    print(f"Reference events:       {reference.size}")
    print(f"Observed events:        {observed.size}")
    print(f"After removefp:         {cleaned.size}")
    print(f"After fillgaps:         {corrected.tn.size}")
    print(f"Unresolved intervals:   {np.count_nonzero(np.isnan(corrected.dtn))}")

    figure, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=False)
    axes[0].eventplot(
        [reference, observed, cleaned, corrected.tn],
        lineoffsets=[4, 3, 2, 1],
        linelengths=0.7,
        colors=["black", "tab:red", "tab:orange", "tab:green"],
    )
    axes[0].set_yticks([1, 2, 3, 4])
    axes[0].set_yticklabels(["Corrected", "Cleaned", "Observed", "Reference"])
    axes[0].set_xlabel("Time (s)")
    axes[0].set_title("Original, detected, cleaned, and reconstructed events")

    axes[1].plot(np.diff(observed), "o-", color="tab:red", label="Observed")
    axes[1].plot(corrected.dtn, "o-", color="tab:green", label="Corrected")
    axes[1].set_xlabel("Interval index")
    axes[1].set_ylabel("Interval (s)")
    axes[1].set_title("Intervals before and after explicit preprocessing")
    axes[1].legend()
    for axis in axes:
        axis.grid(True)
    figure.tight_layout()
    show_or_save_figure(figure, args.save_figure)


if __name__ == "__main__":
    main()
