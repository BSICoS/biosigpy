"""Demonstrate missing-event reconstruction with fillgaps."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.hrv import fillgaps


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
    with_gaps = np.delete(reference, [10, 11, 22])
    corrected = fillgaps(with_gaps)

    print("HRV missing-event reconstruction:")
    print("=================================")
    print(f"Reference events:       {reference.size}")
    print(f"Input events:           {with_gaps.size}")
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
