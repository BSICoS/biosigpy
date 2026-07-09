# `nan_filtfilt`

Import path: `biosigpy.tools.nan_filtfilt.nan_filtfilt`

Apply zero-phase filtering to a one-dimensional signal while preserving long `NaN` gaps.

```python
nan_filtfilt(b, a, x, max_gap=0)
```

## Parameters

| Name | Description |
| --- | --- |
| `b` | Numerator filter coefficients. |
| `a` | Denominator filter coefficients. |
| `x` | One-dimensional real numeric signal. |
| `max_gap` | Maximum internal `NaN` gap length to interpolate before filtering. |

## Returns

A one-dimensional NumPy array with the same length as `x`.

## Notes

Boundary `NaN` gaps are preserved. Internal short `NaN` gaps up to `max_gap` are interpolated for filtering. Internal long `NaN` gaps are preserved and split the signal into separate finite segments. Segments that are too short for zero-phase filtering return `NaN`.

`nan_filtfilt` is conformant with the Biosiglib `tools.nan_filtfilt` specification.

## Example

```python
import numpy as np

from biosigpy.tools import nan_filtfilt

b = np.array([0.25, 0.5, 0.25])
a = np.array([1.0])
x = np.linspace(0.0, 1.0, 20)
x[8] = np.nan

y = nan_filtfilt(b, a, x, max_gap=1)
```
