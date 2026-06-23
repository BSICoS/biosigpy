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
    module, algorithm, case_name = case_id.split(".", maxsplit=2)
    case_path = (
        biosiglib_root()
        / "conformance"
        / module
        / algorithm
        / f"{case_name}.json"
    )
    case_definition = json.loads(case_path.read_text(encoding="utf-8"))
    if case_definition["id"] != case_id:
        raise RuntimeError(f"Conformance case ID mismatch in {case_path}")
    return case_definition


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
