"""Reusable Matplotlib figure handling for examples."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure


_NON_INTERACTIVE_BACKENDS = {
    "agg",
    "cairo",
    "pdf",
    "pgf",
    "ps",
    "svg",
    "template",
}


def add_save_figure_argument(parser: argparse.ArgumentParser) -> None:
    """Add the common ``--save-figure`` option to an example parser."""

    parser.add_argument(
        "--save-figure",
        type=Path,
        default=None,
        metavar="PATH",
        help="Save the generated figure to PATH instead of requiring a GUI.",
    )


def show_or_save_figure(
    figure: Figure,
    save_path: Path | None,
    *,
    dpi: int = 150,
) -> None:
    """Save a figure or show it when an interactive backend is available."""

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(save_path, dpi=dpi)
        print(f"Saved figure: {save_path}")
        return

    backend = plt.get_backend()
    if _is_non_interactive_backend(backend):
        figure.canvas.draw()
        print(
            f"Matplotlib backend {backend!r} is non-interactive; no plot "
            "window was opened. Re-run with --save-figure PATH to save the "
            "figure."
        )
        return

    plt.show()


def _is_non_interactive_backend(backend: str) -> bool:
    normalized_backend = backend.lower()
    if normalized_backend.startswith("module://matplotlib_inline"):
        return True
    return normalized_backend in _NON_INTERACTIVE_BACKENDS
