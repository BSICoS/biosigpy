"""Shared interpolation cases and Python-specific validation."""

import numpy as np
import pytest
from biosigpy.tools import interpgap
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.interpgap"), ids=case_id)
def test_conformance(case_definition):
    result = interpgap(load_input(case_definition, "signal"),
                       load_input(case_definition, "max_gap"),
                       **case_definition["parameters"])
    assert_expected_outputs({"interpolated_signal": result}, case_definition)


@pytest.mark.parametrize("method", ["linear", "nearest", "pchip", "spline"])
@pytest.mark.parametrize("values", [[], [np.nan, np.nan], [1, 2], [np.nan, 2, np.nan]])
def test_identity_paths(values, method):
    np.testing.assert_array_equal(interpgap(values, 2, method), values)


def test_preserves_input_and_normalizes_orientation():
    values = np.array([[0, np.nan, 4]])
    result = interpgap(values, 1)
    np.testing.assert_array_equal(values, [[0, np.nan, 4]])
    np.testing.assert_array_equal(result, [0, 2, 4])


@pytest.mark.parametrize("values, gap, method, error", [
    ([np.inf], 1, "linear", ValueError), ([1j], 1, "linear", TypeError),
    ([[1, 2], [3, 4]], 1, "linear", ValueError),
    ([], -1, "linear", ValueError), ([], 1.5, "linear", TypeError),
    ([], np.inf, "linear", TypeError), ([], True, "linear", TypeError),
    ([], 1, "invalid", ValueError),
])
def test_invalid_inputs(values, gap, method, error):
    with pytest.raises(error):
        interpgap(values, gap, method)
