"""Conformance tests for tools.medfilt_threshold."""

import numpy as np
import pytest

from biosigpy.tools.medfilt_threshold import medfilt_threshold
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    load_case,
    load_input,
)


VALID_CASES = [
    "tools.medfilt_threshold.normal_outlier",
    "tools.medfilt_threshold.max_threshold_cap",
    "tools.medfilt_threshold.window_larger_than_signal",
    "tools.medfilt_threshold.row_vector_orientation",
    "tools.medfilt_threshold.odd_window_behavior",
    "tools.medfilt_threshold.even_window_behavior",
    "tools.medfilt_threshold.include_nan_window",
]

EXPECTED_ERROR_CASES = [
    "tools.medfilt_threshold.invalid_window_one",
    "tools.medfilt_threshold.invalid_single_sample",
]


@pytest.mark.parametrize("case_id", VALID_CASES)
def test_positive_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    parameters = case_definition["parameters"]

    threshold = medfilt_threshold(
        load_input(case_definition, "x"),
        parameters["window"],
        parameters["factor"],
        parameters["max_threshold"],
    )

    assert np.issubdtype(threshold.dtype, np.number)
    assert threshold.ndim == 1
    assert_expected_outputs({"threshold": threshold}, case_definition)


@pytest.mark.parametrize("case_id", EXPECTED_ERROR_CASES)
def test_expected_error_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    parameters = case_definition["parameters"]

    assert_expected_error(
        lambda: medfilt_threshold(
            load_input(case_definition, "x"),
            parameters["window"],
            parameters["factor"],
            parameters["max_threshold"],
        ),
        case_definition,
    )


@pytest.mark.parametrize(
    ("kwargs", "exception_type"),
    [
        ({"x": [[1.0, 2.0], [3.0, 4.0]]}, ValueError),
        ({"x": [1.0, np.inf]}, ValueError),
        ({"window": 2.5}, TypeError),
        ({"factor": 0.0}, ValueError),
        ({"max_threshold": 0.0}, ValueError),
    ],
)
def test_python_validation(
    kwargs: dict[str, object], exception_type: type[Exception]
) -> None:
    arguments = {
        "x": [1.0, 2.0, 3.0],
        "window": 2,
        "factor": 2.0,
        "max_threshold": 10.0,
    }
    arguments.update(kwargs)

    with pytest.raises(exception_type):
        medfilt_threshold(**arguments)
