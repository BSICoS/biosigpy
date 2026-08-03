"""Conformance and Python-native tests for hrv.removefp."""

import numpy as np
import pytest

from biosigpy import removefp as top_level_removefp
from biosigpy.hrv import removefp as hrv_removefp
from biosigpy.hrv.removefp import removefp
from conformance import (
    assert_expected_outputs,
    case_id,
    cases_for_specification,
    load_input,
)


@pytest.mark.parametrize(
    "case_definition",
    cases_for_specification("hrv.removefp"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    result = removefp(load_input(case_definition, "tk"))

    assert result.ndim == 1
    assert result.dtype == np.float64
    assert_expected_outputs({"tn": result}, case_definition)


def test_public_reexports_and_vector_orientation() -> None:
    tk = np.asarray([0, 1, 2, 2.2, 3, 4, 5], dtype=np.float64)
    expected = removefp(tk)

    np.testing.assert_array_equal(hrv_removefp(tk.reshape(1, -1)), expected)
    np.testing.assert_array_equal(top_level_removefp(tk.reshape(-1, 1)), expected)


@pytest.mark.parametrize("tk", ([0, 2, 1], [0, 1, 1], [0, np.nan, 1]))
def test_invalid_event_series_is_rejected(tk: list[float]) -> None:
    with pytest.raises(ValueError):
        removefp(tk)
