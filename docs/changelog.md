# Changelog

## 0.0.0 (in development)

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
- Noted that Biosigpy remains in active development.
