"""Conformance tests for tools.nan_filter."""

import numpy as np
import pytest

from biosigpy.tools.nan_filter import nan_filter
from conformance import assert_expected_outputs, load_case, load_input


VALID_CASES = [
    "tools.nan_filter.no_nan_equivalent_filter",
    "tools.nan_filter.short_nan_gap_interpolation",
    "tools.nan_filter.long_nan_gap_segmentation",
    "tools.nan_filter.row_vector_orientation",
    "tools.nan_filter.boundary_nan_preserved",
    "tools.nan_filter.too_short_segments_nan",
]


@pytest.mark.parametrize("case_id", VALID_CASES)
def test_positive_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    parameters = case_definition["parameters"]

    filtered_signal = nan_filter(
        load_input(case_definition, "numerator_coefficients"),
        load_input(case_definition, "denominator_coefficients"),
        load_input(case_definition, "signal"),
        parameters.get("max_gap", 0),
    )

    assert np.issubdtype(filtered_signal.dtype, np.number)
    assert filtered_signal.ndim == 1
    assert_expected_outputs({"filtered_signal": filtered_signal}, case_definition)


def test_empty_signal_returns_empty() -> None:
    filtered_signal = nan_filter([1.0], [1.0], [])

    assert np.issubdtype(filtered_signal.dtype, np.number)
    assert filtered_signal.ndim == 1
    assert filtered_signal.size == 0


def test_all_nan_signal_returns_all_nan() -> None:
    filtered_signal = nan_filter([1.0], [1.0], [np.nan, np.nan])

    np.testing.assert_equal(filtered_signal, np.array([np.nan, np.nan]))


@pytest.mark.parametrize(
    ("kwargs", "exception_type"),
    [
        ({"b": [1.0, np.inf]}, ValueError),
        ({"a": [1.0, np.inf]}, ValueError),
        ({"x": [1.0, np.inf]}, ValueError),
        ({"x": [[1.0, 2.0], [3.0, 4.0]]}, ValueError),
        ({"max_gap": -1}, ValueError),
        ({"max_gap": 1.5}, TypeError),
    ],
)
def test_python_validation(
    kwargs: dict[str, object], exception_type: type[Exception]
) -> None:
    arguments = {
        "b": [1.0],
        "a": [1.0],
        "x": [1.0, 2.0, 3.0],
        "max_gap": 0,
    }
    arguments.update(kwargs)

    with pytest.raises(exception_type):
        nan_filter(**arguments)
