---
title: ECG Module Overview
---

# ECG Module

The `biosigpy.ecg` package contains conformant ECG preprocessing, detection,
and derived-signal algorithms.

## Functions

- [`baselineremove`](baselineremove.md): estimate and remove ECG baseline
  wander from fiducial isoelectric samples
- [`pantompkins`](pantompkins.md): detect ECG R waves and expose processing
  signals
- [`sloperange`](sloperange.md): derive a respiration signal from ECG slopes

## Examples

- [Baseline-wander removal](../../examples/ecg-baselineremove.md)
- [Pan-Tompkins R-wave detection](../../examples/ecg-pantompkins.md)
- [Slope-range respiration](../../examples/ecg-sloperange.md)

## See also

- [API reference](../index.md)
- [Examples](../../examples/index.md)
