"""Shared tidal-volume cases and Python-specific input validation."""

import numpy as np
import pytest

from biosigpy.resp import tidalvolume
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case", cases_for_specification("resp.tidalvolume"), ids=case_id)
def test_conformance(case):
    resp = load_input(case, "resp")
    fs = load_input(case, "sampling_frequency")
    v, u, l = tidalvolume(resp, fs, **case["parameters"])
    assert v.shape == u.shape == l.shape == (len(resp),)
    assert_expected_outputs(dict(tdvol=v, upper_envelope=u, lower_envelope=l), case)


@pytest.mark.parametrize("signal, fs, kwargs", [
    ([], 1, {}), ([np.inf], 1, {}), ([0], 0, {}),
    ([0], 1, {"mindist": -1}), ([0], 1, {"max_gap": 0}),
])
def test_invalid_inputs(signal, fs, kwargs):
    with pytest.raises(ValueError):
        tidalvolume(signal, fs, **kwargs)
