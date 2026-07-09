# `snap_to_peak`

Import path: `biosigpy.tools.snap_to_peak.snap_to_peak`

Refine approximate ECG detections to local signal maxima.

```python
snap_to_peak(ecg, detections, window_size=20)
```

## Parameters

| Name | Description |
| --- | --- |
| `ecg` | One-dimensional real numeric ECG signal with at least two samples. `NaN` values are allowed and act as finite-segment boundaries. Infinite values are invalid. |
| `detections` | One-dimensional real numeric vector of approximate detection positions. Empty input is allowed. `NaN` detections return `NaN` in the corresponding output position. Infinite values are invalid. |
| `window_size` | Positive search radius in samples. The value is rounded to the nearest integer before constructing search windows. |

## Returns

A one-dimensional NumPy array aligned with `detections`.

## Notes

The public Biosiglib-compatible API uses one-based sample coordinates for `detections` and returned positions. Python internals may use zero-based indices, but callers should pass and interpret this function in one-based sample coordinates.

Search windows are clipped to signal boundaries and to the finite ECG segment containing the detection. If a finite detection lands on a `NaN` ECG sample, the corresponding output is `NaN`. When multiple samples share the maximum value inside the search window, the first maximum is returned.

`snap_to_peak` is conformant with the Biosiglib `tools.snap_to_peak` specification.

## Example

```python
import numpy as np

from biosigpy.tools import snap_to_peak

ecg = np.array([0.1, 0.4, 1.0, 0.5, 0.2])
detections = np.array([2.0])

print(snap_to_peak(ecg, detections, window_size=2))
# [3.]
```
