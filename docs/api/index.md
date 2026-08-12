# API

The API pages document Python signatures and link to executable examples. For expected inputs and outputs, scientific interpretation, limitations, and references, use the [Biosiglib method catalog](https://bsicos.github.io/biosiglib/methods/).

## ECG

| Function | Description |
| --- | --- |
| [`baselineremove`](ecg/baselineremove.md) | Estimate and remove baseline wander from fiducial isoelectric samples. |
| [`pantompkins`](ecg/pantompkins.md) | Detect ECG R waves and expose the processing signals. |
| [`sloperange`](ecg/sloperange.md) | Derive a respiration signal from ECG slopes. |

## HRV

| Function | Description |
| --- | --- |
| [`fdmetrics`](hrv/fdmetrics.md) | Calculate frequency-domain variability metrics. |
| [`fillgaps`](hrv/fillgaps.md) | Reconstruct missing events in an ordered series. |
| [`ipfm`](hrv/ipfm.md) | Reconstruct instantaneous heart rate and optional TVIPFM modulation. |
| [`osp`](hrv/osp.md) | Separate respiration-related HRV modulation from its orthogonal residual. |
| [`removefp`](hrv/removefp.md) | Remove false-positive events from an ordered series. |
| [`tdmetrics`](hrv/tdmetrics.md) | Calculate time-domain variability metrics. |

## Tools

| Function | Description |
| --- | --- |
| [`lpd_filter`](tools/lpd-filter.md) | Apply the low-pass differentiator used by detection algorithms. |
| [`medfilt_threshold`](tools/medfilt-threshold.md) | Compute a median-filtered adaptive threshold. |
| [`nan_filter`](tools/nan-filter.md) | Filter data while preserving NaN gaps. |
| [`nan_filtfilt`](tools/nan-filtfilt.md) | Apply zero-phase filtering while preserving NaN gaps. |
| [`snap_to_peak`](tools/snap-to-peak.md) | Refine event positions by snapping them to nearby peaks. |
