"""Inclusive spectral classification and matching-shape semantics."""

import numpy as np
import pytest

from biosigpy.tools import ispeaky
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.ispeaky"), ids=case_id)
def test_conformance(case_definition):
    pkl = np.asarray(load_input(case_definition, "pkl"), dtype=float)
    akl = np.asarray(load_input(case_definition, "akl"), dtype=float)
    actual = ispeaky(pkl, akl, **case_definition["parameters"])
    assert_expected_outputs(dict(is_peaky=actual.reshape(-1, order="F")), case_definition)


def test_shapes_scalars_and_nan():
    assert ispeaky(45, 85, 45, 85) is True
    assert ispeaky(np.nan, 85, 45, 85) is False
    np.testing.assert_array_equal(ispeaky([[45, 30]], [[85, 100]], 45, 85), [[True, False]])
    for pkl, akl in (([45], [[85]]), ([45, 60], 85), ([], [])):
        with pytest.raises(ValueError):
            ispeaky(pkl, akl, 45, 85)


@pytest.mark.parametrize("threshold", [-1, 101, np.nan, np.inf])
def test_invalid_thresholds(threshold):
    with pytest.raises(ValueError):
        ispeaky([45], [85], threshold, 85)
    with pytest.raises(ValueError):
        ispeaky([45], [85], 45, threshold)
