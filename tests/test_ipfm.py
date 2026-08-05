"""Conformance and Python-native tests for hrv.ipfm."""

import numpy as np
import pytest

from biosigpy import IpfmResult as top_level_result
from biosigpy import ipfm as top_level_ipfm
from biosigpy.hrv import IpfmResult as hrv_result
from biosigpy.hrv import ipfm as hrv_ipfm
from biosigpy.hrv.ipfm import IpfmResult, ipfm
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
    cases_for_specification("hrv.ipfm"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    requested_outputs = case_definition["requested_outputs"]
    return_m = "m" in requested_outputs
    parameters = case_definition["parameters"]

    def execute() -> np.ndarray | IpfmResult:
        return ipfm(
            load_input(case_definition, "tn"),
            load_input(case_definition, "fs"),
            parameters.get("spline_order", 14),
            return_m=return_m,
        )

    if is_expected_error(case_definition):
        assert_expected_error(execute, case_definition)
        return

    result = execute()
    outputs = result._asdict() if isinstance(result, IpfmResult) else {"ihr": result}
    assert_expected_outputs(outputs, case_definition)


def test_named_result_and_public_reexports() -> None:
    result = top_level_ipfm([0, 1, 2, 3], 4, return_m=True)
    ihr, m = result

    assert isinstance(result, top_level_result)
    assert top_level_result is hrv_result is IpfmResult
    assert top_level_ipfm is hrv_ipfm is ipfm
    np.testing.assert_array_equal(ihr, result.ihr)
    np.testing.assert_array_equal(m, result.m)


def test_sampling_frequency_depends_on_requested_outputs() -> None:
    ihr = ipfm([0, 1, 2, 3], 0.05)

    assert np.all(np.isfinite(ihr) & (ihr > 0))
    with pytest.raises(ValueError, match="greater than 0.06"):
        ipfm([0, 1, 2, 3], 0.05, return_m=True)


@pytest.mark.parametrize(
    ("arguments", "keywords", "exception_type"),
    (
        (([], 4), {}, ValueError),
        (([0, 1, 1], 4), {}, ValueError),
        (([0, np.nan], 4), {}, ValueError),
        (([0, 1], "4"), {}, TypeError),
        (([0, 1], 4, 2.0), {}, TypeError),
        (([0, 1], 4), {"return_m": 1}, TypeError),
    ),
)
def test_invalid_inputs_are_rejected(
    arguments: tuple[object, ...],
    keywords: dict[str, object],
    exception_type: type[Exception],
) -> None:
    with pytest.raises(exception_type):
        ipfm(*arguments, **keywords)
