# `pantompkins`

Import path: `biosigpy.ecg.pantompkins.pantompkins`

Detect ECG R waves with a Pan-Tompkins-style processing chain.

```python
pantompkins(
    ecg,
    sampling_frequency,
    *,
    bandpass_frequency=(5.0, 12.0),
    integration_window_size=0.15,
    minimum_peak_distance=0.5,
    snap_to_peak_window_size=20.0,
)
```

## Parameters

| Name | Description |
| --- | --- |
| `ecg` | One-dimensional real numeric ECG signal with at least two samples. |
| `sampling_frequency` | Positive sampling frequency in Hz. |
| `bandpass_frequency` | Two-element bandpass frequency vector in Hz. Values must lie inside the Nyquist range. |
| `integration_window_size` | Moving-integration window length in seconds. |
| `minimum_peak_distance` | Minimum distance between envelope peaks in seconds. |
| `snap_to_peak_window_size` | Search radius in samples used when refining detections to local ECG maxima. |

## Returns

A dictionary with:

| Key | Description |
| --- | --- |
| `r_wave_times` | Detected R-wave times in seconds. |
| `ecg_filtered` | Bandpass-filtered ECG signal. |
| `decg` | Squared derivative signal. |
| `decg_envelope` | Moving-integrated envelope. |

## Notes

`pantompkins` is conformant with the Biosiglib `ecg.pantompkins` specification. It delegates peak refinement to `biosigpy.tools.snap_to_peak.snap_to_peak` and derivative filter design to `biosigpy.tools.lpd_filter.lpd_filter`.

## Example

```python
import numpy as np

from biosigpy.ecg import pantompkins

sampling_frequency = 256.0
time = np.arange(0.0, 2.0, 1.0 / sampling_frequency)
ecg = np.sin(2.0 * np.pi * 1.2 * time)

outputs = pantompkins(ecg, sampling_frequency)
print(outputs["r_wave_times"])
```
