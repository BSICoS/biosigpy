# Changelog

## Unreleased

- Implemented the shared `ecg.baselineremove` contract with structured warning
  metadata, explicit finite-real validation, canonical fiducial normalization,
  local averaging, and MATLAB-compatible spline interpolation and extrapolation.
- Added all eight shared conformance cases, public re-exports, API
  documentation, and a direct Python port of the Biosigmat baseline-removal
  example.
- Added dedicated HRV and ECG module overviews linking every conformant
  function and its executable example, following the Biosigmat structure.
- Implemented the sampled Biosiglib `hrv.ipfm` contract with the canonical
  edge-stabilized B-spline, uniform output grid, and optional TVIPFM signal.
- Consumed all six shared `hrv.ipfm` cases directly from Biosiglib, including
  the Biosigmat-backed nonuniform-event and TVIPFM fixture.
- Added an executable IPFM example aligned with the Biosigmat workflow and
  using the same first 100 Medicom MTD event times.
- Pinned conformance to Biosiglib v1.3.0 at commit
  `a1e7c3eae157c588705117ea9480b2cf3190d107`.

## 0.2.0 - 2026-08-04

- Implemented the Biosiglib `hrv.removefp` contract with strict input order,
  fixed adaptive-baseline settings, and simultaneous one-pass removal.
- Implemented the Biosiglib `hrv.fillgaps` contract with a named two-field
  result, configurable empirical thresholds, segment-wide iterative PCHIP
  reconstruction, exact duration preservation, over-insertion fallback, and
  NaN-marked unresolved gaps.
- Added optional interactive `fillgaps` debugging that displays every
  reconstruction attempt and waits for explicit user inspection without
  changing normative results or affecting normal calls.
- Aligned the Python `removefp` and `fillgaps` examples with the identical
  Medicom MTD timing fixture and transformations used by Biosigmat.
- Added automatic shared-case coverage, Python-native API tests, public
  re-exports, API documentation, and dedicated executable examples for
  `removefp` and `fillgaps`.
- Pinned conformance to Biosiglib v1.2.1 at commit
  `050a0527741415d4099ed6e2d8ba873fc76cf577`.

## 0.1.0 - 2026-08-03

- Published the first reproducible pre-1.0 Biosigpy release as a wheel and
  source distribution with SHA-256 checksums.
- Added one authoritative implementation-version source, installed package
  metadata, and `biosigpy.__version__`.
- Documented and tested Python 3.10-3.13, NumPy 1.26.4+, and SciPy 1.11.4+.
- Added the initial documentation scaffold.
- Documented the current Biosiglib v1.0.0 conformant algorithm set.
- Aligned `hrv.tdmetrics` with the Biosiglib minimum-data contract, returning
  `NaN` only for metrics that are undefined for the available valid intervals.
- Discovered and executed every shared case for each conformant specification,
  with a collection gate that detects unexecuted cases added by Biosiglib.
- Implemented the complete `ecg.sloperange` contract with a named, unpackable
  five-field result, signal-aligned slope vectors, zero-based extrema positions,
  boundary `NaN` values, earliest-sample tie handling, shared conformance cases,
  public API documentation, and an executable visual-inspection example.
- Reached conformance for all eight specifications published by Biosiglib v1.0.0.
- Pinned that contract to Biosiglib v1.0.0 at commit
  `ea2d5ded43e1348342f0db4fbd97d754b90a28c9`.
