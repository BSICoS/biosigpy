"""Conformance and Python-native tests for hrv.fillgaps."""

import numpy as np
import pytest

from biosigpy import FillGapsResult as top_level_result
from biosigpy import fillgaps as top_level_fillgaps
from biosigpy.hrv import FillGapsResult as hrv_result
from biosigpy.hrv import fillgaps as hrv_fillgaps
from biosigpy.hrv.fillgaps import FillGapsResult, fillgaps
from conformance import (
    assert_expected_outputs,
    case_id,
    cases_for_specification,
    load_input,
)


@pytest.mark.parametrize(
    "case_definition",
    cases_for_specification("hrv.fillgaps"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    result = fillgaps(load_input(case_definition, "tk"))

    assert isinstance(result, FillGapsResult)
    assert result.tn.ndim == result.dtn.ndim == 1
    assert result.tn.dtype == result.dtn.dtype == np.float64
    assert_expected_outputs(result._asdict(), case_definition)


def test_named_result_and_public_reexports() -> None:
    result = top_level_fillgaps(np.asarray([[0, 1, 2, 4, 5, 6]]))
    tn, dtn = result

    assert isinstance(result, top_level_result)
    assert top_level_result is hrv_result is FillGapsResult
    assert top_level_fillgaps is hrv_fillgaps is fillgaps
    np.testing.assert_array_equal(tn, result.tn)
    np.testing.assert_array_equal(dtn, result.dtn)


def test_fillgaps_does_not_remove_close_events_implicitly() -> None:
    tk = np.asarray([0, 1, 1.2, 2.2, 3.2, 4.2])

    result = fillgaps(tk)

    assert 1.2 in result.tn
    assert all(event in result.tn for event in tk)


def test_interactive_debug_shows_each_attempt_without_changing_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from biosigpy.hrv import _fillgaps_debug

    pauses: list[int] = []
    monkeypatch.setattr(
        _fillgaps_debug,
        "_is_interactive_backend",
        lambda _backend: True,
    )
    monkeypatch.setattr(
        plt,
        "waitforbuttonpress",
        lambda **_kwargs: pauses.append(1),
    )
    tk = np.asarray(
        [0, 0.4, 1.1, 2.4, 2.9, 3.55], dtype=np.float64
    )
    expected = fillgaps(tk)

    actual = fillgaps(tk, debug=True)

    np.testing.assert_array_equal(actual.tn, expected.tn)
    np.testing.assert_array_equal(actual.dtn, expected.dtn)
    # One rejected single-insertion attempt, one rejected over-insertion,
    # and the retained preceding reconstruction are each inspected.
    assert pauses == [1, 1, 1]
    assert "fillgaps interactive debug" not in plt.get_figlabels()


def test_debug_must_be_boolean() -> None:
    with pytest.raises(TypeError, match="debug must be a boolean"):
        fillgaps([0, 1, 2], debug=1)


def test_debug_requires_an_interactive_backend() -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)

    with pytest.raises(RuntimeError, match="interactive Matplotlib backend"):
        fillgaps([0, 1, 2, 4, 5, 6], debug=True)


@pytest.mark.parametrize(
    "parameters",
    (
        {"gap_detection_factor": np.inf},
        {"correction_lower_factor": 1.15},
        {"correction_upper_factor": 1.6},
        {"minimum_interval": -0.1},
        {"max_gap_duration": 0},
    ),
)
def test_invalid_parameters_are_rejected(parameters: dict[str, float]) -> None:
    with pytest.raises(ValueError):
        fillgaps([0, 1, 2], **parameters)
