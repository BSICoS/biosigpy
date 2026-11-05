"""Shared highest-peak and axis cases with native index origins."""

import numpy as np
import pytest

from biosigpy.tools import localmax
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.localmax"), ids=case_id)
def test_conformance(case_definition):
    signal = load_input(case_definition, "signal")
    canonical_axis = case_definition["parameters"].get("axis")
    axis = None if canonical_axis is None else canonical_axis - 1
    if "expected_error" in case_definition:
        with pytest.raises(ValueError):
            localmax(signal, axis=axis)
        return
    values, positions = localmax(signal, axis=axis)
    assert_expected_outputs(dict(peak_values=np.atleast_1d(values),
                                 peak_positions=np.atleast_1d(positions) + 1), case_definition)


def test_native_axis_and_output_shapes():
    matrix = [[0, 3, 0]]
    values, positions = localmax(matrix)
    np.testing.assert_array_equal(values, [3])
    np.testing.assert_array_equal(positions, [1])
    assert localmax([0, 3, 0]) == (3, 1)
    with pytest.raises(ValueError):
        localmax(matrix, axis=2)
