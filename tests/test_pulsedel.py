"""Shared delineation cases and Python output conventions."""

import warnings

import numpy as np
import pytest

from biosigpy.ppg import pulsedel
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("ppg.pulsedel"), ids=case_id)
def test_conformance(case_definition):
    inputs = {key: load_input(case_definition, key)
              for key in ("ppg", "sampling_frequency", "t_d", "dppg")}
    outputs = pulsedel(**inputs, **case_definition["parameters"])
    assert_expected_outputs(dict(zip(("t_a", "t_b", "t_m", "t_o"), outputs)), case_definition)


@pytest.mark.parametrize("shape", [(3,), (1, 3), (3, 1)])
def test_missing_derivative_and_detection_alignment(shape):
    ppg = [3, 2, 1, 0, 1, 2, 3, 4, 5, 4, 3, 2]
    detections = np.array([.25, np.nan, .25]).reshape(shape)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        outputs = pulsedel(ppg, 20, detections)
    for output, expected in zip(outputs, [.4, .15, .275, np.nan]):
        assert output.shape == (3,)
        np.testing.assert_allclose(output, [expected, np.nan, expected])


def test_empty_detections():
    assert all(x.shape == (0,) for x in pulsedel([1, 2], 20, []))


@pytest.mark.parametrize("time", [-.001, .051, np.inf])
def test_bounds_before_rounding(time):
    with pytest.raises(ValueError, match="within the signal"):
        pulsedel([1, 2], 20, time)

@pytest.mark.parametrize("points", ["A", "b", "M", "o", "BA", "OM", "OAMB"])
def test_selected_points_match_default(points):
    ppg = [3, 2, 1, 0, 1, 2, 3, 4, 5, 4, 3, 2]
    dppg = [0, 0, -1, -.4, .5, 10, 9, 3, -1, -2, -1, 0]
    default = dict(zip("ABMO", pulsedel(ppg, 20, [.25, np.nan], dppg)))
    actual = pulsedel(ppg, 20, [.25, np.nan], dppg, points=points)
    if len(points) == 1:
        assert isinstance(actual, np.ndarray)
        actual = (actual,)
    for value, key in zip(actual, points.upper(), strict=True):
        np.testing.assert_equal(value, default[key])


@pytest.mark.parametrize("points", ["", "AA", "X", None])
def test_invalid_points(points):
    with pytest.raises((TypeError, ValueError)):
        pulsedel([1, 2, 3], 10, [], points=points)
