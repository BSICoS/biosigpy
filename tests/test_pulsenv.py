"""Shared envelope cases and Python input conventions."""

import numpy as np
import pytest

from biosigpy.ppg import pulsedel, pulsenv
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("ppg.pulsenv"), ids=case_id)
def test_conformance(case_definition):
    inputs = {i["id"]: load_input(case_definition, i["id"]) for i in case_definition["inputs"]}
    lower, upper = pulsenv(**inputs, **case_definition["parameters"])
    assert lower.shape == upper.shape == (len(inputs["ppg"]),)
    assert_expected_outputs(dict(lower_envelope=lower, upper_envelope=upper), case_definition)


def test_supplied_fiducials_match_detection_mode():
    x = [3, 2, 1, 0, 1, 2, 3, 4, 5, 4, 3, 2] * 2
    td = np.array([.85, np.nan, .25, .85])
    ta, tb = pulsedel(x, 20, td, points="AB")
    expected = pulsenv(x, 20, td)
    actual = pulsenv(x, 20, t_a=ta, t_b=tb)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)


def test_no_support_and_invalid_inputs():
    with pytest.warns(UserWarning, match="No pulse detections"):
        assert np.isnan(pulsenv([1, 2, 3], 1, [])).all()
    with pytest.warns(UserWarning) as caught:
        assert np.isnan(pulsenv([1, 2, 3], 1, t_a=[1], t_b=[1])).all()
    assert len(caught) == 2
    with pytest.raises(ValueError, match="must be scalar"):
        pulsenv([1, 2, 3], [0, 1, 2], [])
    with pytest.raises(ValueError, match="within the signal"):
        pulsenv([1, 2, 3], 1, [-.01])
    with pytest.raises(ValueError, match="both"):
        pulsenv([1, 2, 3], 1, t_a=[1])
