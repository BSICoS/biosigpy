# Contributing

## Set up a development checkout

Clone Biosigpy next to Biosiglib, then create a virtual environment and install the development dependencies:

```bash
python -m venv .venv
python -m pip install -e ".[dev,docs,examples]"
```

The commit in `biosiglib.lock` is the scientific contract used by the test suite.

## Run the checks

```bash
python -m pytest
```

Every change must preserve the pinned Biosiglib specifications and shared conformance cases. Add or update focused Python tests for implementation-specific behavior.

## Documentation

Public API documentation comes from Python docstrings. Put runnable examples under `examples/` and link them from the corresponding API page.

```bash
python -m mkdocs build --strict
```

## Project rules

- Keep changes small and focused.
- Use English for code, comments, filenames, and technical documentation.
- Preserve scientific formulas, units, defaults, edge cases, and NaN behavior.
- Do not copy specifications, fixtures, or shared cases from Biosiglib.
- Open a pull request only after the complete test suite passes.

Before tagging a release, move its entries from `Unreleased` to a dated version section in `CHANGELOG.md`. The [release workflow](.github/workflows/release.yml) uses that section as the GitHub release notes.
