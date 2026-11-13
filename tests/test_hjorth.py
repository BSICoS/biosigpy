"""Shared Hjorth power, difference-energy and edge cases."""

import pytest

from biosigpy.tools import hjorth
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.hjorth"), ids=case_id)
def test_conformance(case_definition):
    signal = load_input(case_definition, "signal")
    fs = load_input(case_definition, "sampling_frequency")
    if "expected_error" in case_definition:
        with pytest.raises(ValueError):
            hjorth(signal, fs)
        return
    h0, h1, h2 = hjorth(signal, fs)
    assert_expected_outputs(dict(activity=h0, mobility=h1, complexity=h2), case_definition)
