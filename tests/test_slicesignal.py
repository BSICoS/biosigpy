"""Shared window selection, padding and time-center cases."""

import pytest

from biosigpy.tools import slicesignal
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.slicesignal"), ids=case_id)
def test_conformance(case_definition):
    x, window, overlap = (load_input(case_definition, field) for field in
                          ("signal", "window_length", "overlap"))
    params = case_definition["parameters"]
    if "expected_error" in case_definition:
        with pytest.raises(ValueError):
            slicesignal(x, window, overlap, **params)
        return
    result = slicesignal(x, window, overlap, **params)
    actual = {}
    if "sampling_frequency" in params:
        windows, actual["time_centers"] = result
    else:
        windows = result
    assert windows.ndim == 2 and windows.shape[1] == window
    actual["window_samples"] = windows.reshape(-1)
    assert_expected_outputs(actual, case_definition)
