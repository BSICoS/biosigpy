# Changelog

This file is the authoritative release history. GitHub Releases reuse the corresponding version section.

## 0.3.0 - 2026-08-12

- Simplified the public documentation to Home, Installation, and API; moved contributor guidance to `CONTRIBUTING.md` and linked examples directly from function pages.
- Replaced hardcoded package versions in installation documentation with the latest GitHub release.
- Replaced `conformance.json` with a one-line `biosiglib.lock` against one exact Biosiglib commit; removed per-specification statuses, entry points, and duplicated release and repository metadata.
- Discover every specification and shared case in the pinned Biosiglib commit, rejecting specifications without cases and uncollected cases without a support filter.
- Changed the Biosiglib pin workflow to validate a commit before release rather than propagating partial support after release.
- Implemented the shared `ecg.baselineremove` contract with structured warning metadata, explicit finite-real validation, canonical fiducial normalization, local averaging, and MATLAB-compatible spline interpolation and extrapolation.
- Added all eight shared conformance cases, public re-exports, API documentation, and a direct Python port of the Biosigmat baseline-removal example.
- Implemented the sampled Biosiglib `hrv.ipfm` contract with the canonical edge-stabilized B-spline, uniform output grid, and optional TVIPFM signal.
- Consumed all six shared `hrv.ipfm` cases directly from Biosiglib, including the Biosigmat-backed nonuniform-event and TVIPFM fixture.
- Added an executable IPFM example aligned with the Biosigmat workflow and using the same first 100 Medicom MTD event times.
- Pinned conformance to the Biosiglib v2.0.0 release target at commit `9f2370451c6b77296c2714f719132edb7fa034c0`.

## 0.2.0 - 2026-08-04

[GitHub release](https://github.com/BSICoS/biosigpy/releases/tag/v0.2.0)

- Implemented the Biosiglib `hrv.removefp` contract with strict input order, fixed adaptive-baseline settings, and simultaneous one-pass removal.
- Implemented the Biosiglib `hrv.fillgaps` contract with a named two-field result, configurable empirical thresholds, segment-wide iterative PCHIP reconstruction, exact duration preservation, over-insertion fallback, and NaN-marked unresolved gaps.
- Added optional interactive `fillgaps` debugging that displays every reconstruction attempt and waits for explicit user inspection without changing normative results or affecting normal calls.
- Aligned the Python `removefp` and `fillgaps` examples with the identical Medicom MTD timing fixture and transformations used by Biosigmat.
- Added automatic shared-case coverage, Python-native API tests, public re-exports, API documentation, and dedicated executable examples for `removefp` and `fillgaps`.
- Pinned conformance to Biosiglib v1.2.1 at commit `050a0527741415d4099ed6e2d8ba873fc76cf577`.

## 0.1.0 - 2026-08-03

[GitHub release](https://github.com/BSICoS/biosigpy/releases/tag/v0.1.0)

- Published the first reproducible pre-1.0 Biosigpy release as a wheel and source distribution with SHA-256 checksums.
- Added one authoritative implementation-version source, installed package metadata, and `biosigpy.__version__`.
- Documented and tested Python 3.10–3.13, NumPy 1.26.4+, and SciPy 1.11.4+.
- Added the initial documentation scaffold.
- Documented the current Biosiglib v1.0.0 conformant algorithm set.
- Aligned `hrv.tdmetrics` with the Biosiglib minimum-data contract, returning `NaN` only for metrics that are undefined for the available valid intervals.
- Discovered and executed every shared case for each conformant specification, with a collection gate that detects unexecuted cases added by Biosiglib.
- Implemented the complete `ecg.sloperange` contract with a named, unpackable five-field result, signal-aligned slope vectors, zero-based extrema positions, boundary `NaN` values, earliest-sample tie handling, shared conformance cases, public API documentation, and an executable visual-inspection example.
- Reached conformance for all eight specifications published by Biosiglib v1.0.0.
- Pinned that contract to Biosiglib v1.0.0 at commit `ea2d5ded43e1348342f0db4fbd97d754b90a28c9`.
