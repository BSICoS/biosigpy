"""Conformance and Python-native tests for hrv.fillgaps."""

import numpy as np
import pytest

from biosigpy import FillGapsResult as top_level_result
from biosigpy import fillgaps as top_level_fillgaps
from biosigpy.hrv import FillGapsResult as hrv_result
from biosigpy.hrv import fillgaps as hrv_fillgaps
from biosigpy.hrv.fillgaps import FillGapsResult, fillgaps
from conformance import (
    assert_expected_outputs,
    case_id,
    cases_for_specification,
    load_input,
)


@pytest.mark.parametrize(
    "case_definition",
    cases_for_specification("hrv.fillgaps"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    result = fillgaps(load_input(case_definition, "tk"))

    assert isinstance(result, FillGapsResult)
    assert result.tn.ndim == result.dtn.ndim == 1
    assert result.tn.dtype == result.dtn.dtype == np.float64
    assert_expected_outputs(result._asdict(), case_definition)


def test_named_result_and_public_reexports() -> None:
    result = top_level_fillgaps(np.asarray([[0, 1, 2, 4, 5, 6]]))
    tn, dtn = result

    assert isinstance(result, top_level_result)
    assert top_level_result is hrv_result is FillGapsResult
    assert top_level_fillgaps is hrv_fillgaps is fillgaps
    np.testing.assert_array_equal(tn, result.tn)
    np.testing.assert_array_equal(dtn, result.dtn)


def test_fillgaps_does_not_remove_close_events_implicitly() -> None:
    tk = np.asarray([0, 1, 1.2, 2.2, 3.2, 4.2])

    result = fillgaps(tk)

    assert 1.2 in result.tn
    assert all(event in result.tn for event in tk)


@pytest.mark.parametrize(
    "parameters",
    (
        {"gap_detection_factor": np.inf},
        {"correction_lower_factor": 1.15},
        {"correction_upper_factor": 1.6},
        {"minimum_interval": -0.1},
        {"max_gap_duration": 0},
    ),
)
def test_invalid_parameters_are_rejected(parameters: dict[str, float]) -> None:
    with pytest.raises(ValueError):
        fillgaps([0, 1, 2], **parameters)

