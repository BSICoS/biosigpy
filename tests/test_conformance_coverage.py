"""Regression tests for shared conformance harness behavior."""

from types import SimpleNamespace

import pytest

import conftest
from conformance import assert_complete_case_coverage, assert_expected_warnings


def test_newly_discovered_uncollected_case_fails() -> None:
    discovered = {
        "hrv.tdmetrics.existing_case",
        "hrv.tdmetrics.new_case",
    }
    collected = {"hrv.tdmetrics.existing_case"}

    with pytest.raises(RuntimeError, match="hrv.tdmetrics.new_case"):
        assert_complete_case_coverage(discovered, collected)


@pytest.mark.parametrize("raw_args", ((), ("-q",)))
def test_full_suite_still_rejects_an_uncollected_case(
    raw_args: tuple[str, ...],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    existing_case = {
        "id": "hrv.tdmetrics.existing_case",
        "specification_id": "hrv.tdmetrics",
        "expected_error": {"category": "invalid_value"},
    }
    missing_case = {
        "id": "hrv.tdmetrics.new_case",
        "specification_id": "hrv.tdmetrics",
        "expected_error": {"category": "invalid_value"},
    }
    monkeypatch.setattr(
        conftest,
        "discover_cases",
        lambda: [existing_case, missing_case],
    )
    config = _fake_config(raw_args=raw_args)
    item = SimpleNamespace(
        callspec=SimpleNamespace(params={"case_definition": existing_case})
    )

    with pytest.raises(pytest.UsageError, match="hrv.tdmetrics.new_case"):
        conftest.pytest_collection_modifyitems(config, [item])


@pytest.mark.parametrize(
    ("args_source", "option_overrides"),
    (
        ("ARGS", {}),
        ("INVOCATION_DIR", {"keyword": "selected"}),
        ("INVOCATION_DIR", {"markexpr": "slow"}),
        ("INVOCATION_DIR", {"deselect": ["tests/test_example.py"]}),
        ("INVOCATION_DIR", {"lf": True}),
    ),
)
def test_scoped_collections_skip_the_global_gate(
    args_source: str,
    option_overrides: dict[str, object],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        conftest,
        "discover_cases",
        lambda: [{"id": "hrv.tdmetrics.uncollected"}],
    )
    config = _fake_config(
        args_source=args_source,
        option_overrides=option_overrides,
    )

    conftest.pytest_collection_modifyitems(config, [])


def test_expected_warnings_compare_ids_and_affected_sets() -> None:
    case_definition = {
        "expected_warnings": [
            {"id": "example_warning", "affected_ids": ["first", "second"]}
        ]
    }
    actual = [
        SimpleNamespace(
            warning_id="example_warning",
            affected_ids=("second", "first"),
        )
    ]

    assert_expected_warnings(actual, case_definition)


def test_duplicate_canonical_warning_id_fails() -> None:
    case_definition = {
        "expected_warnings": [
            {"id": "example_warning", "affected_ids": ["output"]}
        ]
    }
    warning = SimpleNamespace(
        warning_id="example_warning",
        affected_ids=("output",),
    )

    with pytest.raises(AssertionError, match="once per canonical id"):
        assert_expected_warnings([warning, warning], case_definition)


def _fake_config(
    *,
    raw_args: tuple[str, ...] = (),
    args_source: str = "INVOCATION_DIR",
    option_overrides: dict[str, object] | None = None,
) -> SimpleNamespace:
    options: dict[str, object] = {
        "keyword": "",
        "markexpr": "",
        "ignore": None,
        "ignore_glob": None,
        "deselect": None,
        "lf": False,
    }
    options.update(option_overrides or {})
    return SimpleNamespace(
        args_source=SimpleNamespace(name=args_source),
        invocation_params=SimpleNamespace(args=raw_args),
        option=SimpleNamespace(**options),
    )
