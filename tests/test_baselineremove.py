"""Conformance and Python-native tests for ecg.baselineremove."""

from __future__ import annotations

import warnings
from collections.abc import Callable

import numpy as np
import pytest

from biosigpy import baselineremove as top_level_baselineremove
from biosigpy.ecg import baselineremove as ecg_baselineremove
from biosigpy.ecg.baselineremove import (
    BaselineRemoveResult,
    BaselineRemoveWarning,
    baselineremove,
)
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    assert_expected_warnings,
    case_id,
    cases_for_specification,
    is_expected_error,
    load_input,
)


@pytest.mark.parametrize(
    "case_definition",
    cases_for_specification("ecg.baselineremove"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    run_case = lambda: _execute_case(case_definition)
    if is_expected_error(case_definition):
        assert_expected_error(run_case, case_definition)
        return

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_case()

    actual_warnings = [
        item.message
        for item in caught
        if isinstance(item.message, BaselineRemoveWarning)
    ]
    assert_expected_warnings(actual_warnings, case_definition)
    assert len(actual_warnings) == len(caught)
    assert result.ecg_detrended.ndim == 1
    assert result.baseline.ndim == 1
    assert_expected_outputs(result._asdict(), case_definition)


def _execute_case(
    case_definition: dict[str, object],
) -> BaselineRemoveResult:
    parameters = case_definition["parameters"]
    assert isinstance(parameters, dict)
    return baselineremove(
        load_input(case_definition, "ecg"),
        load_input(case_definition, "fiducial_positions"),
        load_input(case_definition, "offset"),
        parameters.get("window_size", 5),
    )


def test_row_and_column_vectors_have_identical_semantics() -> None:
    ecg = np.array([2.0, 4.0, 8.0, 16.0, 32.0])
    positions = np.array([1.0, 5.0])

    expected = baselineremove(ecg, positions, 0, 2)
    actual = baselineremove(
        ecg.reshape(1, -1), positions.reshape(-1, 1), 0, 2
    )

    np.testing.assert_array_equal(actual.ecg_detrended, expected.ecg_detrended)
    np.testing.assert_array_equal(actual.baseline, expected.baseline)


@pytest.mark.parametrize(
    "execute",
    (
        lambda: baselineremove([1 + 1j, 2], [1, 2], 0),
        lambda: baselineremove([1, 2], [1 + 1j, 2], 0),
        lambda: baselineremove([1, 2], [1, 2], -1),
        lambda: baselineremove([1, 2], [1, 2], 0, 0),
    ),
)
def test_python_validation(execute: Callable[[], object]) -> None:
    with pytest.raises((TypeError, ValueError)):
        execute()


def test_named_result_warning_and_public_reexports() -> None:
    with pytest.warns(BaselineRemoveWarning) as caught:
        result = baselineremove([1, 2, 3], [10], 0)

    assert isinstance(result, BaselineRemoveResult)
    assert result._fields == ("ecg_detrended", "baseline")
    assert caught[0].message.warning_id == "no_valid_fiducial_positions"
    assert caught[0].message.affected_ids == ("fiducial_positions",)
    assert ecg_baselineremove is baselineremove
    assert top_level_baselineremove is baselineremove
