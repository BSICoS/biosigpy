"""Conformance tests for tools.medfilt_threshold."""

import numpy as np
import pytest

from biosigpy.tools.medfilt_threshold import medfilt_threshold
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
    cases_for_specification("tools.medfilt_threshold"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    parameters = case_definition["parameters"]
    run_case = lambda: medfilt_threshold(
        load_input(case_definition, "x"),
        parameters["window"],
        parameters["factor"],
        parameters["max_threshold"],
    )
    if is_expected_error(case_definition):
        assert_expected_error(run_case, case_definition)
        return

    threshold = run_case()

    assert np.issubdtype(threshold.dtype, np.number)
    assert threshold.ndim == 1
    assert_expected_outputs({"threshold": threshold}, case_definition)


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
