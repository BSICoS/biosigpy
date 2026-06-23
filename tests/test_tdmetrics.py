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


EXPECTED_ERROR_CASES = [
    "hrv.tdmetrics.invalid_tk_non_numeric",
    "hrv.tdmetrics.invalid_tk_matrix",
    "hrv.tdmetrics.invalid_tk_non_monotonic",
    "hrv.tdmetrics.invalid_tk_repeated",
    "hrv.tdmetrics.invalid_tk_negative",
]


def test_positive_conformance() -> None:
    case_definition = load_case("hrv.tdmetrics.ecg_tk_001")
    outputs = tdmetrics(load_input(case_definition, "tk"))

    assert set(outputs) == {"mhr", "sdnn", "sdsd", "rmssd", "pnn50"}
    assert_expected_outputs(outputs, case_definition)


@pytest.mark.parametrize("case_id", EXPECTED_ERROR_CASES)
def test_expected_error_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    tk = load_input(case_definition, "tk")

    assert_expected_error(lambda: tdmetrics(tk), case_definition)


def test_row_and_column_vectors_match_one_dimensional_input() -> None:
    tk = np.array([0.0, 0.8, 1.7, 2.5, 3.4])
    expected = tdmetrics(tk)

    assert tdmetrics(tk.reshape(1, -1)) == expected
    assert tdmetrics(tk.reshape(-1, 1)) == expected
