"""Shared conformance and Python input tests for tools.trimnans."""

import numpy as np
import pytest

from biosigpy.tools import trimnans
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.trimnans"), ids=case_id)
def test_conformance(case_definition):
    result = trimnans(load_input(case_definition, "signal"))
    assert result.ndim == 1
    assert_expected_outputs({"trimmed_signal": result}, case_definition)


@pytest.mark.parametrize("signal", ([], [np.nan, np.nan]))
def test_no_retained_samples(signal):
    result = trimnans(signal)
    assert result.ndim == 1
    assert result.size == 0


@pytest.mark.parametrize("signal, error", [
    ([np.inf], ValueError), ([-np.inf], ValueError),
    ([np.nan, np.inf, 2], ValueError), ([1+2j, 3], TypeError),
    ("invalid", TypeError), ([[1, 2], [3, 4]], ValueError),
])
def test_invalid_inputs(signal, error):
    with pytest.raises(error):
        trimnans(signal)


@pytest.mark.parametrize("shape", ((1, 4), (4, 1)))
def test_vector_orientation(shape):
    signal = np.array([np.nan, 1, 2, np.nan]).reshape(shape)
    np.testing.assert_array_equal(trimnans(signal), [1, 2])
