"""Helpers for consuming pinned Biosiglib conformance resources."""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from collections.abc import Callable, Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any, Protocol

import numpy as np
import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class CanonicalWarning(Protocol):
    """Warning data required by the shared conformance comparison."""

    warning_id: str
    affected_ids: Sequence[str]


@lru_cache(maxsize=1)
def load_biosiglib_commit() -> str:
    """Read the exact Biosiglib commit pinned by this repository."""

    lock_path = REPOSITORY_ROOT / "biosiglib.lock"
    try:
        raw_lock = lock_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise RuntimeError(f"Unable to read Biosiglib lock {lock_path}: {exc}") from exc
    if re.fullmatch(r"[0-9a-f]{40}\n?", raw_lock) is None:
        raise RuntimeError(
            f"Biosiglib lock must contain one lowercase 40-character SHA: {lock_path}"
        )
    return raw_lock.strip()


@lru_cache(maxsize=1)
def biosiglib_root() -> Path:
    configured_root = os.environ.get("BIOSIGLIB_ROOT")
    root = (
        Path(configured_root)
        if configured_root
        else REPOSITORY_ROOT.parent / "biosiglib"
    ).resolve()
    if not root.is_dir():
        raise RuntimeError(f"Biosiglib checkout does not exist: {root}")

    expected_commit = load_biosiglib_commit()
    completed = subprocess.run(
        [
            "git",
            "-c",
            f"safe.directory={root.as_posix()}",
            "-C",
            str(root),
            "rev-parse",
            "HEAD",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    actual_commit = completed.stdout.strip()
    if actual_commit != expected_commit:
        raise RuntimeError(
            f"Biosiglib commit mismatch: expected {expected_commit}, got {actual_commit}"
        )
    return root


@lru_cache(maxsize=None)
def load_case(case_id: str) -> dict[str, Any]:
    matches = [case for case in discover_cases() if case["id"] == case_id]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one discovered conformance case {case_id!r}; "
            f"found {len(matches)}"
        )
    return matches[0]


@lru_cache(maxsize=1)
def discover_cases() -> tuple[dict[str, Any], ...]:
    """Discover every case for every specification in the pinned checkout."""

    discovered: list[dict[str, Any]] = []
    seen_case_ids: set[str] = set()

    for specification_id in discover_specification_ids():
        module, algorithm = specification_id.split(".", maxsplit=1)
        case_directory = biosiglib_root() / "conformance" / module / algorithm
        case_paths = sorted(case_directory.glob("*.json"))
        if not case_paths:
            raise RuntimeError(
                "No Biosiglib conformance cases found for "
                f"specification {specification_id!r} in {case_directory}"
            )

        for case_path in case_paths:
            case_definition = json.loads(case_path.read_text(encoding="utf-8"))
            case_id = case_definition.get("id")
            if case_definition.get("specification_id") != specification_id:
                raise RuntimeError(
                    f"Conformance case specification mismatch in {case_path}"
                )
            expected_case_id = f"{specification_id}.{case_path.stem}"
            if case_id != expected_case_id:
                raise RuntimeError(
                    f"Conformance case ID mismatch in {case_path}: "
                    f"expected {expected_case_id!r}, got {case_id!r}"
                )
            if case_id in seen_case_ids:
                raise RuntimeError(f"Duplicate conformance case ID {case_id!r}")

            has_outputs = "expected_outputs" in case_definition
            has_error = "expected_error" in case_definition
            if has_outputs == has_error:
                raise RuntimeError(
                    f"Conformance case {case_id!r} must define exactly one of "
                    "expected_outputs or expected_error"
                )
            seen_case_ids.add(case_id)
            discovered.append(case_definition)

    return tuple(discovered)


@lru_cache(maxsize=1)
def discover_specification_ids() -> tuple[str, ...]:
    """Return every canonical specification ID from the pinned checkout."""

    specification_ids: set[str] = set()
    for specification_path in sorted((biosiglib_root() / "specs").rglob("spec.json")):
        specification = json.loads(specification_path.read_text(encoding="utf-8"))
        specification_id = specification.get("metadata", {}).get("id")
        if not isinstance(specification_id, str) or not specification_id:
            raise RuntimeError(
                f"Biosiglib specification has no canonical ID: {specification_path}"
            )
        if specification_id in specification_ids:
            raise RuntimeError(
                f"Duplicate Biosiglib specification ID {specification_id!r}"
            )
        specification_ids.add(specification_id)

    if not specification_ids:
        raise RuntimeError("No Biosiglib specifications found in the pinned checkout")
    return tuple(sorted(specification_ids))


def cases_for_specification(specification_id: str) -> tuple[dict[str, Any], ...]:
    """Return all discovered cases for one pinned specification."""

    if specification_id not in discover_specification_ids():
        raise RuntimeError(f"Unknown Biosiglib specification {specification_id!r}")
    return tuple(
        case
        for case in discover_cases()
        if case["specification_id"] == specification_id
    )


def case_id(case_definition: Mapping[str, Any]) -> str:
    """Return the canonical ID used for a parametrized-test label."""

    return str(case_definition["id"])


def is_expected_error(case_definition: Mapping[str, Any]) -> bool:
    """Classify a case from its Biosiglib definition."""

    return "expected_error" in case_definition


def missing_case_ids(
    discovered_case_ids: set[str], collected_case_ids: set[str]
) -> set[str]:
    """Return discovered cases that are absent from the test collection."""

    return discovered_case_ids - collected_case_ids


def assert_complete_case_coverage(
    discovered_case_ids: set[str], collected_case_ids: set[str]
) -> None:
    """Reject a test collection that omits any discovered case."""

    missing = missing_case_ids(discovered_case_ids, collected_case_ids)
    if missing:
        formatted = ", ".join(sorted(missing))
        raise RuntimeError(
            "Discovered Biosiglib conformance cases were not collected: "
            + formatted
        )


@lru_cache(maxsize=1)
def load_fixture_catalog() -> dict[str, Any]:
    catalog_path = biosiglib_root() / "fixtures" / "catalog.json"
    return json.loads(catalog_path.read_text(encoding="utf-8"))


@lru_cache(maxsize=None)
def load_fixture_column(fixture_id: str, file_role: str, column: str) -> np.ndarray:
    fixture = _exactly_one(
        load_fixture_catalog()["fixtures"], "id", fixture_id, "fixture"
    )
    fixture_file = _exactly_one(
        fixture["files"], "role", file_role, "fixture file role"
    )
    if fixture_file["format"] != "csv":
        raise RuntimeError(f"Fixture {fixture_id} role {file_role} is not CSV")

    root = biosiglib_root()
    csv_path = (root / fixture_file["path"]).resolve()
    if not csv_path.is_relative_to(root):
        raise RuntimeError(f"Fixture path escapes Biosiglib: {csv_path}")

    with csv_path.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise RuntimeError(f"Fixture column {column!r} does not exist in {csv_path}")
        return np.asarray([float(row[column]) for row in reader], dtype=np.float64)


def load_input(case_definition: Mapping[str, Any], input_id: str) -> Any:
    input_definition = _exactly_one(
        case_definition["inputs"], "id", input_id, "case input"
    )
    if "value" in input_definition:
        return input_definition["value"]
    return load_fixture_column(
        input_definition["fixture_id"],
        input_definition["file_role"],
        input_definition["column"],
    )


def assert_expected_outputs(
    actual_outputs: Mapping[str, Any], case_definition: Mapping[str, Any]
) -> None:
    for expected_output in case_definition["expected_outputs"]:
        output_id = expected_output["id"]
        assert output_id in actual_outputs, f"Missing canonical output {output_id!r}"
        expected_value = (
            expected_output["value"]
            if "value" in expected_output
            else load_fixture_column(
                expected_output["fixture_id"],
                expected_output["file_role"],
                expected_output["column"],
            )
        )
        expected_value = _decode_expected_value(expected_value)
        np.testing.assert_allclose(
            actual_outputs[output_id],
            expected_value,
            rtol=0.0,
            atol=expected_output["absolute_tolerance"],
            equal_nan=case_definition["nan_equal"],
        )


def assert_expected_warnings(
    actual_warnings: Sequence[CanonicalWarning],
    case_definition: Mapping[str, Any],
) -> None:
    """Compare canonical warning IDs and aggregated affected IDs."""

    expected_warnings = case_definition.get("expected_warnings", [])
    actual_by_id = {warning.warning_id: warning for warning in actual_warnings}
    expected_by_id = {
        expected_warning["id"]: expected_warning
        for expected_warning in expected_warnings
    }

    assert len(actual_warnings) == len(actual_by_id), (
        "warnings must be emitted once per canonical id"
    )
    assert actual_by_id.keys() == expected_by_id.keys()
    for warning_id, expected_warning in expected_by_id.items():
        assert set(actual_by_id[warning_id].affected_ids) == set(
            expected_warning["affected_ids"]
        )


def _decode_expected_value(value: Any) -> Any:
    if isinstance(value, str):
        if value != "NaN":
            raise RuntimeError(f"Unsupported expected string value {value!r}")
        return np.nan
    if isinstance(value, list):
        return [_decode_expected_value(item) for item in value]
    return value


def assert_expected_error(
    function: Callable[[], Any], case_definition: Mapping[str, Any]
) -> None:
    category = case_definition["expected_error"]["category"]
    exception_types: dict[str, type[Exception] | tuple[type[Exception], ...]] = {
        "invalid_type": (TypeError, ValueError),
        "invalid_shape": ValueError,
        "invalid_value": ValueError,
        "insufficient_data": ValueError,
        "invalid_numerical_result": ValueError,
    }
    with pytest.raises(exception_types[category]):
        function()


def _exactly_one(
    items: list[dict[str, Any]], field: str, value: str, description: str
) -> dict[str, Any]:
    matches = [item for item in items if item.get(field) == value]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one {description} with {field}={value!r}; "
            f"found {len(matches)}"
        )
    return matches[0]
