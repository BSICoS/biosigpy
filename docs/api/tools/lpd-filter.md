# `lpd_filter`

Import path: `biosigpy.tools.lpd_filter.lpd_filter`

Design a MATLAB/Biosigmat-compatible low-pass differentiator FIR filter.

```python
lpd_filter(sampling_frequency, stop_frequency, pass_frequency=None, order=None)
```

## Parameters

| Name | Description |
| --- | --- |
| `sampling_frequency` | Positive sampling frequency in Hz. |
| `stop_frequency` | Positive stop frequency in Hz, below Nyquist. |
| `pass_frequency` | Optional positive pass frequency in Hz. When omitted, Biosigpy uses `stop_frequency - 0.2`. |
| `order` | Explicit positive filter order. Automatic order selection is currently unsupported. |

## Returns

A tuple `(b, delay)` where `b` is a one-dimensional NumPy array of FIR coefficients and `delay` is the filter delay in samples.

## Notes

`lpd_filter` is conformant with the Biosiglib `tools.lpd_filter` specification. It supports explicit order design only. If `order` is `None`, the function raises `NotImplementedError`.

The requested order is rounded up to an even value. The returned delay is `filter_order / 2`, where `filter_order` is the rounded filter order used for coefficient design.

## Example

```python
from biosigpy.tools import lpd_filter

b, delay = lpd_filter(256.0, 12.0, order=4)
```
