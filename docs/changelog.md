# Changelog

## 0.0.0 (in development)

- Added the initial documentation scaffold.
- Documented the current Biosiglib v0.6.0 conformant algorithm set.
- Aligned `hrv.tdmetrics` with the Biosiglib minimum-data contract, returning
  `NaN` only for metrics that are undefined for the available valid intervals.
- Discovered and executed every shared case for each conformant specification,
  with a collection gate that detects unexecuted cases added by Biosiglib.
- Implemented `ecg.sloperange`, including aligned boundary `NaN` values,
  MATLAB-compatible sample-grid rounding, shared conformance cases, public API
  documentation, and an executable ECG-derived respiration example.
- Reached conformance for all eight specifications published by Biosiglib v0.6.0.
- Noted that Biosigpy remains in active development.
