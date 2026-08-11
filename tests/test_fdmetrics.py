"""Conformance and Python-native tests for hrv.fdmetrics."""

from __future__ import annotations

import warnings
from collections.abc import Callable

import numpy as np
import pytest

from biosigpy import fdmetrics as top_level_fdmetrics
from biosigpy.hrv import fdmetrics as hrv_fdmetrics
from biosigpy.hrv.fdmetrics import (
    FdMetricsResult,
    FdMetricsWarning,
    SeparatedFdMetricsResult,
    fdmetrics,
)
from conformance import (
    assert_expected_outputs,
    assert_expected_warnings,
    case_id,
    cases_for_specification,
    load_input,
)


@pytest.mark.parametrize(
    "case_definition",
    cases_for_specification("hrv.fdmetrics"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = _execute_case(case_definition)

    actual_warnings = [
        item.message
        for item in caught
        if isinstance(item.message, FdMetricsWarning)
    ]
    assert_expected_warnings(actual_warnings, case_definition)
    assert len(actual_warnings) == len(caught)
    assert_expected_outputs(result._asdict(), case_definition)


def _execute_case(
    case_definition: dict[str, object],
) -> FdMetricsResult | SeparatedFdMetricsResult:
    input_ids = {item["id"] for item in case_definition["inputs"]}
    if "pxx" in input_ids:
        return fdmetrics(
            load_input(case_definition, "pxx"),
            load_input(case_definition, "f"),
            case_definition["parameters"].get("limit_hf", True),
        )
    return fdmetrics(
        f=load_input(case_definition, "f"),
        related_pxx=load_input(case_definition, "related_pxx"),
        unrelated_pxx=load_input(case_definition, "unrelated_pxx"),
    )


def test_named_results_and_public_reexports() -> None:
    conventional = fdmetrics([1, 1, 1], [0.04, 0.15, 0.4])
    separated = fdmetrics(
        f=[0.04, 0.15, 0.4],
        related_pxx=[0.01, 0.01, 0.01],
        unrelated_pxx=[0.001, 0.001, 0.001],
    )

    assert isinstance(conventional, FdMetricsResult)
    assert isinstance(separated, SeparatedFdMetricsResult)
    assert top_level_fdmetrics is hrv_fdmetrics is fdmetrics


@pytest.mark.parametrize(
    "execute",
    (
        lambda: fdmetrics([], []),
        lambda: fdmetrics([1, -1], [0.04, 0.15]),
        lambda: fdmetrics([1, 1], [0.04, np.inf]),
        lambda: fdmetrics([1, 1], [0.15, 0.04]),
        lambda: fdmetrics([1], [0.04, 0.15]),
        lambda: fdmetrics([1, 1], [0.04, 0.15], 1),
        lambda: fdmetrics(
            [1, 1],
            [0.04, 0.15],
            related_pxx=[1, 1],
            unrelated_pxx=[1, 1],
        ),
        lambda: fdmetrics(f=[0.04, 0.15], related_pxx=[1, 1]),
    ),
)
def test_invalid_inputs_are_rejected(execute: Callable[[], object]) -> None:
    with pytest.raises((TypeError, ValueError)):
        execute()
