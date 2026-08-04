"""Interactive visual inspection for :mod:`biosigpy.hrv.fillgaps`."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray


class FillGapsDebugger:
    """Render each fill attempt and wait for explicit user interaction."""

    _figure_label = "fillgaps interactive debug"

    def __init__(self) -> None:
        try:
            import matplotlib.pyplot as plt
        except ImportError as error:
            raise ImportError(
                "debug=True requires Matplotlib; install biosigpy[examples]"
            ) from error
        if not _is_interactive_backend(plt.get_backend()):
            raise RuntimeError(
                "debug=True requires an interactive Matplotlib backend"
            )

        self._plt = plt
        self._figure = None
        self._overview_axis = None
        self._attempt_axis = None

    def overview(
        self,
        intervals: NDArray[np.float64],
        gap_indices: NDArray[np.int64],
        detection_threshold: NDArray[np.float64],
        insertion_count: int,
    ) -> None:
        """Show the current intervals, correctable gaps, and threshold."""

        self._ensure_figure()
        axis = self._overview_axis
        assert axis is not None

        axis.clear()
        indices = np.arange(intervals.size)
        _stem(
            axis,
            indices,
            intervals,
            color="0.15",
            label="RR intervals",
        )
        _stem(
            axis,
            gap_indices,
            intervals[gap_indices],
            color="tab:red",
            linewidth=1.5,
            label="Detected gaps",
        )
        axis.plot(
            indices,
            detection_threshold,
            "k--",
            label="Detection threshold",
        )
        axis.set_xlim(-0.5, max(0.5, intervals.size - 0.5))
        axis.set_ylabel("Current RR [s]")
        axis.set_title(
            "Detected gaps: testing "
            f"{insertion_count} inserted event(s)"
        )
        axis.legend(loc="best")
        self._draw()

    def attempt(
        self,
        intervals: NDArray[np.float64],
        gap_index: int,
        insertion_count: int,
        upper_boundary: float,
        lower_boundary: float,
        *,
        accepted: bool,
    ) -> None:
        """Show one candidate, then wait for a key press or mouse click."""

        self._ensure_figure()
        axis = self._attempt_axis
        assert axis is not None

        axis.clear()
        indices = np.arange(intervals.size)
        _stem(
            axis,
            indices,
            intervals,
            color="0.15",
            label="Candidate intervals",
        )
        filled_indices = np.arange(
            gap_index, gap_index + insertion_count + 1
        )
        status = "accepted" if accepted else "rejected"
        color = "tab:green" if accepted else "tab:red"
        _stem(
            axis,
            filled_indices,
            intervals[filled_indices],
            color=color,
            linewidth=1.5,
            label=status,
        )

        left = max(-0.5, gap_index - 50.0)
        right = min(intervals.size - 0.5, gap_index + 50.0)
        axis.set_xlim(left, max(left + 1.0, right))
        visible_maximum = max(
            float(np.max(intervals[filled_indices])),
            upper_boundary,
            lower_boundary,
        )
        if not np.isfinite(visible_maximum) or visible_maximum <= 0:
            visible_maximum = 1.0
        axis.set_ylim(0.0, 1.1 * visible_maximum)
        axis.axhline(
            upper_boundary,
            color="black",
            linestyle="--",
            label="Upper limit",
        )
        axis.axhline(
            lower_boundary,
            color="black",
            linestyle=":",
            label="Lower limit",
        )
        axis.set_xlabel("Interval index (zero-based)")
        axis.set_ylabel("Candidate RR [s]")
        axis.set_title(
            f"{status.capitalize()} attempt: "
            f"{insertion_count} inserted event(s)"
        )
        axis.legend(loc="best")
        self._draw()
        self._plt.waitforbuttonpress(timeout=-1)

    def close(self) -> None:
        """Close the debug figure if an attempt created it."""

        if self._figure is not None:
            self._plt.close(self._figure)

    def _ensure_figure(self) -> None:
        if self._figure is not None:
            return
        self._figure = self._plt.figure(
            num=self._figure_label, figsize=(14, 8)
        )
        self._overview_axis, self._attempt_axis = self._figure.subplots(2, 1)
        self._figure.tight_layout()

    def _draw(self) -> None:
        assert self._figure is not None
        self._figure.tight_layout()
        self._figure.canvas.draw_idle()
        self._figure.show(warn=False)


def _is_interactive_backend(backend: str) -> bool:
    normalized = backend.casefold()
    return normalized not in {
        "agg",
        "cairo",
        "pdf",
        "pgf",
        "ps",
        "svg",
        "template",
    } and not normalized.startswith("module://matplotlib_inline")


def _stem(
    axis: Any,
    indices: NDArray[np.int64],
    values: NDArray[np.float64],
    *,
    color: str,
    linewidth: float = 1.0,
    label: str,
) -> None:
    """Draw a small stem plot without depending on StemContainer styling."""

    axis.vlines(
        indices,
        0.0,
        values,
        colors=color,
        linewidth=linewidth,
        label=label,
    )
    axis.plot(indices, values, "o", color=color)
