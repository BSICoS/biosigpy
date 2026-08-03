"""Helpers for consuming pinned Biosiglib conformance resources."""

from __future__ import annotations

import csv
import json
import os
import subprocess
from collections.abc import Callable, Mapping
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


@lru_cache(maxsize=1)
def load_manifest() -> dict[str, Any]:
    return json.loads((REPOSITORY_ROOT / "conformance.json").read_text(encoding="utf-8"))


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

    expected_commit = load_manifest()["biosiglib"]["commit"]
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
    """Discover every case for specifications declared conformant."""

    discovered: list[dict[str, Any]] = []
    seen_case_ids: set[str] = set()
    manifest = load_manifest()
    conformant_specification_ids = sorted(
        specification_id
        for specification_id, implementation in manifest["specifications"].items()
        if implementation["status"] == "conformant"
    )

    for specification_id in conformant_specification_ids:
        module, algorithm = specification_id.split(".", maxsplit=1)
        case_directory = biosiglib_root() / "conformance" / module / algorithm
        case_paths = sorted(case_directory.glob("*.json"))
        if not case_paths:
            raise RuntimeError(
                "No Biosiglib conformance cases found for conformant "
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


def cases_for_specification(specification_id: str) -> tuple[dict[str, Any], ...]:
    """Return all discovered cases for one conformant specification."""

    manifest_entry = load_manifest()["specifications"].get(specification_id)
    if manifest_entry is None or manifest_entry["status"] != "conformant":
        raise RuntimeError(
            f"Specification {specification_id!r} is not declared conformant"
        )
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
        if isinstance(expected_value, str) and expected_value == "NaN":
            expected_value = np.nan
        np.testing.assert_allclose(
            actual_outputs[output_id],
            expected_value,
            rtol=0.0,
            atol=expected_output["absolute_tolerance"],
            equal_nan=case_definition["nan_equal"],
        )


def assert_expected_error(
    function: Callable[[], Any], case_definition: Mapping[str, Any]
) -> None:
    category = case_definition["expected_error"]["category"]
    exception_types: dict[str, type[Exception] | tuple[type[Exception], ...]] = {
        "invalid_type": (TypeError, ValueError),
        "invalid_shape": ValueError,
        "invalid_value": ValueError,
        "insufficient_data": ValueError,
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
