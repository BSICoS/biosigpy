# API

The API pages document Python signatures and link to executable examples. For expected inputs and outputs, scientific interpretation, limitations, and references, use the [Biosiglib method catalog](https://bsicos.github.io/biosiglib/methods/).

## ECG

| Function | Description |
| --- | --- |
| [`baselineremove`](ecg/baselineremove.md) | Estimate and remove baseline wander from fiducial isoelectric samples. |
| [`pantompkins`](ecg/pantompkins.md) | Detect ECG R waves and expose the processing signals. |
| [`sloperange`](ecg/sloperange.md) | Derive a respiration signal from ECG slopes. |

## PPG

| Function | Description |
| --- | --- |
| [`hjorth_artifacts`](ppg/hjorth-artifacts.md) | Detect artifacts from windowed Hjorth descriptors. |

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
| [`lpdfilter`](tools/lpdfilter.md) | Apply the low-pass differentiator used by detection algorithms. |
| [`medfilt_threshold`](tools/medfilt-threshold.md) | Compute a median-filtered adaptive threshold. |
| [`nanfilter`](tools/nanfilter.md) | Filter data while preserving NaN gaps. |
| [`nanfiltfilt`](tools/nanfiltfilt.md) | Apply zero-phase filtering while preserving NaN gaps. |
| [`snaptopeak`](tools/snaptopeak.md) | Refine event positions by snapping them to nearby peaks. |
| [`trimnans`](tools/trimnans.md) | Trim boundary NaNs while preserving internal gaps. |
| [`expandnans`](tools/expandnans.md) | Expand original NaN gaps by a time radius. |
| [`findsequences`](tools/findsequences.md) | Locate repeated values with zero-based inclusive run boundaries. |
| [`slicesignal`](tools/slicesignal.md) | Slice signals into overlapping windows. |
| [`localmax`](tools/localmax.md) | Find the highest local peak along a signal axis. |
| [`hjorth`](tools/hjorth.md) | Compute mean-square activity, mobility, and frequency spread. |
| [`nanpwelch`](tools/nanpwelch.md) | Estimate PSD from valid stretches without crossing gaps. |
| [`peakedness`](tools/peakedness.md) | Measure spectral peakedness. |
| [`ispeaky`](tools/ispeaky.md) | Classify spectral peakedness. |
