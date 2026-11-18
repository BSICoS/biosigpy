"""Detect PPG artifacts using the same signal and parameters as Biosigmat.

Run from the repository root with:
    python examples/ppg/hjorth_artifacts_example.py
"""

from pathlib import Path

import numpy as np

from biosigpy.ppg import hjorth_artifacts


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"


def main():
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


if __name__ == "__main__":
    main()
