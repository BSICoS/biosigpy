"""Conformance tests for tools.nan_filtfilt."""

import numpy as np
import pytest

from biosigpy.tools.nan_filtfilt import nan_filtfilt
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
    cases_for_specification("tools.nan_filtfilt"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    parameters = case_definition["parameters"]
    run_case = lambda: nan_filtfilt(
        load_input(case_definition, "numerator_coefficients"),
        load_input(case_definition, "denominator_coefficients"),
        load_input(case_definition, "signal"),
        parameters.get("max_gap", 0),
    )
    if is_expected_error(case_definition):
        assert_expected_error(run_case, case_definition)
        return

    filtered_signal = run_case()

    assert np.issubdtype(filtered_signal.dtype, np.number)
    assert filtered_signal.ndim == 1
    assert_expected_outputs({"filtered_signal": filtered_signal}, case_definition)


def test_empty_signal_returns_empty() -> None:
    filtered_signal = nan_filtfilt([1.0], [1.0], [])

    assert np.issubdtype(filtered_signal.dtype, np.number)
    assert filtered_signal.ndim == 1
    assert filtered_signal.size == 0


def test_all_nan_signal_returns_all_nan() -> None:
    filtered_signal = nan_filtfilt([1.0], [1.0], [np.nan, np.nan])

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
        nan_filtfilt(**arguments)
