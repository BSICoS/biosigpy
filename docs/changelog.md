# Changelog

## Unreleased

- No changes yet.

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
