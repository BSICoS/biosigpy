"""Shared run detection with native-to-canonical index normalization."""

import numpy as np
import pytest
from biosigpy.tools import findsequences
from conformance import assert_expected_outputs, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition", cases_for_specification("tools.findsequences"), ids=case_id)
def test_conformance(case_definition):
    result = findsequences(load_input(case_definition, "signal"))
    assert result.ndim == 2 and result.shape[1] == 4
    assert_expected_outputs(dict(values=result[:, 0], start_positions=result[:, 1] + 1,
                                 end_positions=result[:, 2] + 1, lengths=result[:, 3]), case_definition)


def test_native_indices_and_slicing():
    signal = np.array([5, 5, 9, 2, 2, 2])
    result = findsequences(signal)
    np.testing.assert_array_equal(result, [[5, 0, 1, 2], [2, 3, 5, 3]])
    for value, first, last, length in result:
        run = signal[int(first):int(last)+1]
        assert run.size == length
        assert np.all(run == value)


def test_integer_input_and_orientation():
    np.testing.assert_array_equal(findsequences(np.array([2, 1, 1, 0, 0], dtype=np.uint8)),
                                  [[1, 1, 2, 2], [0, 3, 4, 2]])
    np.testing.assert_array_equal(findsequences([[1], [1], [2]]), [[1, 0, 1, 2]])


@pytest.mark.parametrize("signal, error", [([1j, 1j], TypeError), ("text", TypeError),
    ([True, True], TypeError), ([[1, 2], [3, 4]], ValueError)])
def test_invalid_input(signal, error):
    with pytest.raises(error):
        findsequences(signal)
