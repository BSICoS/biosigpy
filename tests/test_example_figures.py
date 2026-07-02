"""Tests for reusable example figure helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from examples._support.figures import (
    add_save_figure_argument,
    show_or_save_figure,
)


def test_add_save_figure_argument_parses_path() -> None:
    parser = argparse.ArgumentParser()
    add_save_figure_argument(parser)

    args = parser.parse_args(["--save-figure", "figure.png"])

    assert args.save_figure == Path("figure.png")


def test_show_or_save_figure_writes_file(tmp_path, capsys) -> None:
    output_path = tmp_path / "nested" / "figure.png"
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])

    show_or_save_figure(figure, output_path)

    assert output_path.is_file()
    assert f"Saved figure: {output_path}" in capsys.readouterr().out
    plt.close(figure)


def test_non_interactive_backend_does_not_call_show(monkeypatch, capsys) -> None:
    def fail_show() -> None:
        raise AssertionError("plt.show() must not be called for Agg")

    monkeypatch.setattr(plt, "show", fail_show)
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])

    show_or_save_figure(figure, None)

    output = capsys.readouterr().out
    assert "non-interactive" in output
    assert "--save-figure PATH" in output
    plt.close(figure)
