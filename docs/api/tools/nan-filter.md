# `nan_filter`

Import path: `biosigpy.tools.nan_filter.nan_filter`

Apply causal filtering to a one-dimensional signal while preserving long `NaN` gaps.

```python
nan_filter(b, a, x, max_gap=0)
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

Boundary `NaN` gaps are preserved. Internal short `NaN` gaps up to `max_gap` are interpolated for filtering. Internal long `NaN` gaps are preserved and split the signal into separate finite segments. Segments that are too short for causal filtering return `NaN`.

`nan_filter` is conformant with the Biosiglib `tools.nan_filter` specification.

## Example

```python
import numpy as np

from biosigpy.tools import nan_filter

b = np.array([0.5, 0.5])
a = np.array([1.0])
x = np.array([1.0, np.nan, 3.0, 4.0])

y = nan_filter(b, a, x, max_gap=1)
```
