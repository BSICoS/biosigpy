"""Shared expansion cases and Python-specific validation."""

import numpy as np
import pytest

from biosigpy.tools import expandnans
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.expandnans"), ids=case_id)
def test_conformance(case_definition):
    result = expandnans(*(load_input(case_definition, field) for field in
                         ("signal", "sampling_frequency", "expansion_duration")))
    assert result.ndim == 1
    assert_expected_outputs({"expanded_signal": result}, case_definition)


@pytest.mark.parametrize("signal, duration", [([], 0.2), ([np.nan, np.nan], 0.2),
    ([1, 2, 3], 0.2), ([1, np.nan, 3], 0), ([1, np.nan, 3], 0.04)])
def test_identity_paths(signal, duration):
    np.testing.assert_array_equal(expandnans(signal, 10, duration), signal)


@pytest.mark.parametrize("signal, fs, duration, error", [
    ([np.inf], 10, 0.2, ValueError), ([-np.inf], 10, 0.2, ValueError),
    ([1+2j], 10, 0.2, TypeError), ([[1, 2], [3, 4]], 10, 0.2, ValueError),
    ([], 0, 0.2, ValueError), ([], np.inf, 0.2, ValueError),
    ([], 10, np.nan, ValueError), ([], 10, -1, ValueError),
    ([], 10, 1j, TypeError), ([], 10, "invalid", TypeError),
])
def test_invalid_inputs(signal, fs, duration, error):
    with pytest.raises(error):
        expandnans(signal, fs, duration)
