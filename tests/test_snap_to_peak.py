"""Conformance tests for tools.snap_to_peak."""

import numpy as np
import pytest

from biosigpy.tools.snap_to_peak import snap_to_peak
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    load_case,
    load_input,
)


VALID_CASES = [
    "tools.snap_to_peak.local_maxima",
    "tools.snap_to_peak.boundary_clipping",
    "tools.snap_to_peak.configurable_window_small",
    "tools.snap_to_peak.configurable_window_large",
    "tools.snap_to_peak.ecg_nan_segment_boundary",
    "tools.snap_to_peak.detection_nan_returns_nan",
    "tools.snap_to_peak.detection_on_nan_ecg_returns_nan",
]

EXPECTED_ERROR_CASES = [
    "tools.snap_to_peak.invalid_detection_out_of_bounds",
]


@pytest.mark.parametrize("case_id", VALID_CASES)
def test_positive_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    parameters = case_definition["parameters"]

    refined_detections = snap_to_peak(
        load_input(case_definition, "ecg"),
        load_input(case_definition, "detections"),
        parameters.get("window_size", 20),
    )

    assert np.issubdtype(refined_detections.dtype, np.number)
    assert refined_detections.ndim == 1
    assert_expected_outputs(
        {"refined_detections": refined_detections}, case_definition
    )


@pytest.mark.parametrize("case_id", EXPECTED_ERROR_CASES)
def test_expected_error_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    parameters = case_definition["parameters"]

    assert_expected_error(
        lambda: snap_to_peak(
            load_input(case_definition, "ecg"),
            load_input(case_definition, "detections"),
            parameters.get("window_size", 20),
        ),
        case_definition,
    )


@pytest.mark.parametrize(
    ("kwargs", "exception_type"),
    [
        ({"ecg": [1.0, np.inf]}, ValueError),
        ({"detections": [np.inf]}, ValueError),
        ({"ecg": [[1.0, 2.0], [3.0, 4.0]]}, ValueError),
    ],
)
def test_python_validation(
    kwargs: dict[str, object], exception_type: type[Exception]
) -> None:
    arguments = {
        "ecg": [1.0, 3.0, 1.0],
        "detections": [2.0],
        "window_size": 1.0,
    }
    arguments.update(kwargs)

    with pytest.raises(exception_type):
        snap_to_peak(**arguments)


def test_empty_detections_returns_empty() -> None:
    refined_detections = snap_to_peak([1.0, 2.0], [])

    assert np.issubdtype(refined_detections.dtype, np.number)
    assert refined_detections.ndim == 1
    assert refined_detections.size == 0


def test_first_maximum_wins_for_ties() -> None:
    refined_detections = snap_to_peak([0.0, 5.0, 5.0, 1.0], [3.0], window_size=2)

    np.testing.assert_array_equal(refined_detections, np.array([2.0]))
