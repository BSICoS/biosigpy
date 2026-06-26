"""Conformance tests for hrv.tdmetrics."""

import numpy as np
import pytest

from biosigpy.hrv.tdmetrics import tdmetrics
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    load_case,
    load_input,
)


VALID_CASES = [
    "hrv.tdmetrics.valid_dtk_001",
    "hrv.tdmetrics.valid_dtk_with_nan_001",
]

EXPECTED_ERROR_CASES = [
    "hrv.tdmetrics.invalid_dtk_non_numeric",
    "hrv.tdmetrics.invalid_dtk_matrix",
    "hrv.tdmetrics.invalid_dtk_negative",
    "hrv.tdmetrics.invalid_dtk_zero",
    "hrv.tdmetrics.invalid_dtk_inf",
]


@pytest.mark.parametrize("case_id", VALID_CASES)
def test_positive_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    outputs = tdmetrics(load_input(case_definition, "dtk"))

    assert set(outputs) == {"mhr", "sdnn", "sdsd", "rmssd", "pnn50"}
    assert_expected_outputs(outputs, case_definition)


@pytest.mark.parametrize("case_id", EXPECTED_ERROR_CASES)
def test_expected_error_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    dtk = load_input(case_definition, "dtk")

    assert_expected_error(lambda: tdmetrics(dtk), case_definition)


def test_nan_markers_are_omitted_before_successive_differences() -> None:
    dtk = np.array([0.8, np.nan, 0.82, 0.78])
    expected = tdmetrics(np.array([0.8, 0.82, 0.78]))

    assert tdmetrics(dtk) == expected
