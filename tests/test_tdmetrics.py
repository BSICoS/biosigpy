"""Conformance tests for hrv.tdmetrics."""

import numpy as np
import pytest

from biosigpy.hrv.tdmetrics import tdmetrics
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
    cases_for_specification("hrv.tdmetrics"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    dtk = load_input(case_definition, "dtk")
    if is_expected_error(case_definition):
        assert_expected_error(lambda: tdmetrics(dtk), case_definition)
        return

    outputs = tdmetrics(dtk)

    assert set(outputs) == {"mhr", "sdnn", "sdsd", "rmssd", "pnn50"}
    assert_expected_outputs(outputs, case_definition)


def test_nan_markers_are_omitted_before_successive_differences() -> None:
    dtk = np.array([0.8, np.nan, 0.82, 0.78])
    expected = tdmetrics(np.array([0.8, 0.82, 0.78]))

    assert tdmetrics(dtk) == expected


def test_all_nan_input_returns_all_nan_metrics() -> None:
    outputs = tdmetrics(np.array([np.nan, np.nan]))

    assert set(outputs) == {"mhr", "sdnn", "sdsd", "rmssd", "pnn50"}
    assert all(np.isnan(value) for value in outputs.values())
