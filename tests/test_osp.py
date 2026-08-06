"""Conformance and Python-native tests for hrv.osp."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest

from biosigpy import OspResult as top_level_result
from biosigpy import osp as top_level_osp
from biosigpy.hrv import OspResult as hrv_result
from biosigpy.hrv import osp as hrv_osp
from biosigpy.hrv.osp import OspResult, osp
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
    cases_for_specification("hrv.osp"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    parameters = case_definition["parameters"]

    def execute() -> OspResult:
        return osp(
            load_input(case_definition, "m"),
            load_input(case_definition, "resp"),
            load_input(case_definition, "resp_pxx"),
            load_input(case_definition, "f"),
            load_input(case_definition, "fs"),
            parameters.get("min_resp_frequency", 0.1),
        )

    if is_expected_error(case_definition):
        assert_expected_error(execute, case_definition)
        return

    assert_expected_outputs(execute()._asdict(), case_definition)


def test_named_result_and_public_reexports() -> None:
    result = osp(
        [99, 1, 2, 3, 4, 5],
        [1, 0, -1, 0, 1, 0],
        [0, 0, 1],
        [0, 0.5, 1],
        1,
    )
    m_resp, m_unrelated, delay = result

    assert isinstance(result, top_level_result)
    assert top_level_result is hrv_result is OspResult
    assert top_level_osp is hrv_osp is osp
    np.testing.assert_array_equal(m_resp, result.m_resp)
    np.testing.assert_array_equal(m_unrelated, result.m_unrelated)
    assert delay == result.delay


@pytest.mark.parametrize(
    "execute",
    (
        lambda: osp([1], [1], [1], [0], 1),
        lambda: osp([1], [1], [1, -1], [0, 1], 1),
        lambda: osp([1], [1], [1, 1], [0, np.nan], 1),
        lambda: osp([1], [1], [1, 1], [1, 0], 1),
        lambda: osp([1], [1], [1, 1], [0, 1], 0),
        lambda: osp([1], [1], [1, 1], [0, 1], 1, 0),
        lambda: osp([1, 2], [1], [1, 1], [0, 1], 1),
        lambda: osp("m", [1], [1, 1], [0, 1], 1),
    ),
)
def test_invalid_inputs_are_rejected(execute: Callable[[], object]) -> None:
    with pytest.raises((TypeError, ValueError)):
        execute()
