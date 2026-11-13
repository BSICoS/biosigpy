"""Shared spectral-band calculations and local API validation."""

import numpy as np
import pytest

from biosigpy.tools import peakedness
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.peakedness"), ids=case_id)
def test_conformance(case_definition):
    pxx = np.asarray(load_input(case_definition, "pxx"), dtype=float)
    f = load_input(case_definition, "f")
    if "expected_error" in case_definition:
        with pytest.raises(ValueError):
            peakedness(pxx, f, **case_definition["parameters"])
        return
    pkl, akl = peakedness(pxx, f, **case_definition["parameters"])
    assert_expected_outputs(dict(pkl=np.atleast_1d(pkl), akl=np.atleast_1d(akl)), case_definition)


@pytest.mark.parametrize("overrides", [
    {"pxx": [1, np.nan, 2]}, {"pxx": [1, np.inf, 2]},
    {"pxx": [1, -1, 2]}, {"pxx": [1, 1j, 2]}, {"pxx": []},
    {"f": [0, 1, 1]}, {"f": [0, 2, 1]}, {"f": [0, np.nan, 2]},
    {"f": [0, np.inf, 2]}, {"f": [0, 1]},
    {"window": 0}, {"window": np.inf}, {"reference_freq": np.nan},
])
def test_invalid_inputs(overrides):
    args = dict(pxx=[1, 3, 1], f=[0, 1, 2])
    args.update(overrides)
    with pytest.raises((ValueError, TypeError)):
        peakedness(**args)


def test_vector_orientation_and_integer_input():
    expected = peakedness([1, 3, 1], [0, 1, 2], None, 4)
    assert peakedness([[1, 3, 1]], [[0, 1, 2]], None, 4) == expected
    assert peakedness([[1], [3], [1]], [0, 1, 2], None, 4) == expected
    assert peakedness(np.array([1, 3, 1], dtype=np.uint8), [0, 1, 2], None, 4) == expected
