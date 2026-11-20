"""Shared pulse-detection cases and Python input conventions."""

import numpy as np
import pytest

from biosigpy.ppg import pulsedet
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize(
    "case_definition", cases_for_specification("ppg.pulsedet"), ids=case_id
)
def test_conformance(case_definition):
    signal = load_input(case_definition, "dppg")
    fs = load_input(case_definition, "sampling_frequency")
    t_d, threshold = pulsedet(signal, fs, **case_definition["parameters"])
    assert t_d.ndim == threshold.ndim == 1
    assert threshold.size == signal.size
    assert_expected_outputs(dict(t_d=t_d, threshold=threshold), case_definition)


@pytest.mark.parametrize("shape", [(5,), (1, 5), (5, 1)])
def test_vector_orientation_and_input_preservation(shape):
    signal = np.array([0, 0, 1, 0, 0]).reshape(shape)
    original = signal.copy()
    t_d, threshold = pulsedet(signal, 10)
    np.testing.assert_allclose(t_d, [.2], atol=1e-14)
    np.testing.assert_allclose(threshold, [.6, .54, 1, 1, 1], atol=1e-14)
    np.testing.assert_array_equal(signal, original)
    assert threshold.shape == (5,)


def test_infinite_reset_disables_optional_feature():
    signal = np.array([0, 0, 1, 0, 0])
    expected = pulsedet(signal, 10)
    actual = pulsedet(signal, 10, reset_after=np.inf)
    for output, reference in zip(actual, expected):
        np.testing.assert_array_equal(output, reference)


def test_rounded_zero_decay_keeps_undefined_endpoint():
    t_d, threshold = pulsedet([0, 0, 1, 0, 0], 10, adaptive_tau_rr=.001)
    np.testing.assert_allclose(t_d, [.2], atol=1e-14)
    np.testing.assert_allclose(threshold, [.6, .54, 1, 1, np.nan], atol=1e-14)


@pytest.mark.parametrize("overrides", [
    {"dppg": []}, {"dppg": [0, np.inf, 0]}, {"dppg": [0, 1j, 0]},
    {"method": "unknown"}, {"reset_after": 0}, {"adaptive_alpha_amp": 0},
])
def test_invalid_inputs(overrides):
    args = dict(dppg=[0, 0, 1, 0, 0], sampling_frequency=10)
    args.update(overrides)
    with pytest.raises((TypeError, ValueError)):
        pulsedet(**args)
