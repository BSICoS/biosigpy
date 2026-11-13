"""Detect PPG artifacts using the same signal and parameters as Biosigmat.

Run from the repository root with:
    python examples/ppg/hjorth_artifacts_example.py
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

from biosigpy.ppg import hjorth_artifacts


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
if str(EXAMPLES_ROOT) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_ROOT))

from _support.figures import add_save_figure_argument, show_or_save_figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_save_figure_argument(parser)
    args = parser.parse_args()

    # Load PPG signal from fixtures
    ppg_data = np.genfromtxt(
        EXAMPLES_ROOT / "fixtures/ppg/ppg_signals.csv", delimiter=",", names=True,
    )
    ppg = ppg_data["sig"]
    fs = 1000

    # Normalize -> This helps to establish margins (MATLAB normalize defaults)
    ppg = (ppg - np.mean(ppg)) / np.std(ppg, ddof=1)

    # Define parameters
    seg = 4
    step = 3
    margin_h0 = [5, 1]
    margin_h1 = [0.5, 0.5]
    margin_h2 = [1, 2]
    margins = [margin_h0, margin_h1, margin_h2]

    # Get both artifact vector and matrix
    artifact_vector, artifact_matrix = hjorth_artifacts(ppg, fs, seg, step, margins, plotflag=True)

    show_or_save_figure(plt.gcf(), args.save_figure)


if __name__ == "__main__":
    main()
