"""Shared window, threshold and interval-merging cases."""

import numpy as np
import pytest
from biosigpy.ppg import hjorth_artifacts
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("ppg.hjorth_artifacts"), ids=case_id)
def test_conformance(case_definition):
    args = [load_input(case_definition, key) for key in
            ("signal","sampling_frequency","window_duration","step_duration","margins")]
    mask, intervals = hjorth_artifacts(*args, **case_definition["parameters"])
    assert mask.dtype == np.bool_
    assert intervals.ndim == 2 and intervals.shape[1] == 2
    assert_expected_outputs(dict(artifact_mask=mask, artifact_starts=intervals[:,0],
                                 artifact_ends=intervals[:,1]), case_definition)



@pytest.mark.parametrize("overrides", [
    {"window_duration": 0.02}, {"step_duration": 0.001},
    {"signal": [0, 1]}, {"min_segment_separation": -1},
    {"median_filter_order": 1.5}, {"margins": [[-1, 1]] * 3},
])
def test_invalid_window_and_options(overrides):
    args = dict(signal=[0, 1, 0] * 4, sampling_frequency=100,
                window_duration=0.03, step_duration=0.03,
                margins=[[1, 1]] * 3)
    args.update(overrides)
    with pytest.raises((TypeError, ValueError)):
        hjorth_artifacts(**args)



@pytest.mark.parametrize("signal", [[0, 1, 0] * 4, [2] * 12, [np.nan] * 12])
def test_plotflag_preserves_results_and_has_four_panels(signal):
    import matplotlib.pyplot as plt

    args = (signal, 1, 3, 3, [[100, 100]] * 3)
    expected_mask, expected_intervals = hjorth_artifacts(*args)
    mask, intervals = hjorth_artifacts(*args, plotflag=True)
    np.testing.assert_array_equal(mask, expected_mask)
    np.testing.assert_array_equal(intervals, expected_intervals)
    figure = plt.gcf()
    try:
        assert [axis.get_ylabel() for axis in figure.axes] == [
            "Input", "Activity", "Mobility", "Complexity"
        ]
        for axis in figure.axes[1:]:
            assert len(axis.lines) == 3  # Descriptor, upper and lower thresholds
    finally:
        plt.close(figure)
