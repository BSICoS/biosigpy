# `medfilt_threshold`

Import path: `biosigpy.tools.medfilt_threshold.medfilt_threshold`

Compute a capped median-filtered adaptive threshold.

```python
medfilt_threshold(x, window, factor, max_threshold)
```

## Parameters

| Name | Description |
| --- | --- |
| `x` | One-dimensional real numeric signal with at least two samples. Infinite values are invalid. |
| `window` | Median window length. Values below 2 are invalid. Windows longer than the signal are clipped to the signal length. |
| `factor` | Positive multiplier applied to the median-filtered signal. |
| `max_threshold` | Positive maximum threshold value. |

## Returns

A one-dimensional NumPy array with the same length as `x`.

## Notes

`medfilt_threshold` is conformant with the Biosiglib `tools.medfilt_threshold` specification.

## Example

```python
import numpy as np

from biosigpy.tools import medfilt_threshold

x = np.array([1.0, 1.2, 2.4, 1.1, 1.0])
threshold = medfilt_threshold(x, window=3, factor=1.5, max_threshold=2.0)
```
