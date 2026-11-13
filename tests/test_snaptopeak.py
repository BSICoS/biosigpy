"""Conformance tests for tools.snaptopeak."""

import numpy as np
import pytest

from biosigpy.tools.snaptopeak import snaptopeak
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    case_id,
    cases_for_specification,
    is_expected_error,
    load_input,
)


@pytest.mark.parametrize(
    "case_definition",
    cases_for_specification("tools.snaptopeak"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    parameters = case_definition["parameters"]
    run_case = lambda: snaptopeak(
        load_input(case_definition, "ecg"),
        np.asarray(load_input(case_definition, "detections"), dtype=float) - 1,
        parameters.get("window_size", 20),
    )
    if is_expected_error(case_definition):
        assert_expected_error(run_case, case_definition)
        return

    refined_detections = run_case()

    assert np.issubdtype(refined_detections.dtype, np.number)
    assert refined_detections.ndim == 1
    assert_expected_outputs(
        {"refined_detections": refined_detections + 1}, case_definition
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
        "detections": [1.0],
        "window_size": 1.0,
    }
    arguments.update(kwargs)

    with pytest.raises(exception_type):
        snaptopeak(**arguments)


def test_empty_detections_returns_empty() -> None:
    refined_detections = snaptopeak([1.0, 2.0], [])

    assert np.issubdtype(refined_detections.dtype, np.number)
    assert refined_detections.ndim == 1
    assert refined_detections.size == 0


def test_first_maximum_wins_for_ties() -> None:
    refined_detections = snaptopeak([0.0, 5.0, 5.0, 1.0], [2.0], window_size=2)

    np.testing.assert_array_equal(refined_detections, np.array([1.0]))


def test_native_boundary_indices():
    result = snaptopeak([1, 0, 0, 2], [0, 3], window_size=0.1)
    np.testing.assert_array_equal(result, [0, 3])


@pytest.mark.parametrize("position", [-0.1, 3.1])
def test_reject_positions_outside_native_grid(position):
    with pytest.raises(ValueError):
        snaptopeak([1, 0, 0, 2], [position])


def test_fractional_positions_and_radius_round_half_up():
    # A radius below half a sample leaves the rounded detection in place.
    result = snaptopeak([4, 2, 3, 5], [0.49, 0.5, 2.5], window_size=0.1)
    np.testing.assert_array_equal(result, [0, 1, 3])
    # A half-sample radius includes the next sample in the search.
    np.testing.assert_array_equal(snaptopeak([0, 4, 0], [0], 0.5), [1])
