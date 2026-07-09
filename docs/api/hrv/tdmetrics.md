# `tdmetrics`

Import path: `biosigpy.hrv.tdmetrics.tdmetrics`

Compute time-domain beat or pulse variability metrics from an interval series.

```python
tdmetrics(dtk)
```

## Parameters

| Name | Description |
| --- | --- |
| `dtk` | One-dimensional real numeric interval series in seconds. `NaN` values are allowed as missing interval markers and are omitted before computing metrics. Infinite, zero, and negative values are invalid. |

## Returns

A `dict[str, float]` with:

| Key | Unit | Description |
| --- | --- | --- |
| `mhr` | beats/min | Mean heart rate. |
| `sdnn` | ms | Standard deviation of intervals. |
| `sdsd` | ms | Standard deviation of successive interval differences. |
| `rmssd` | ms | Root mean square of successive interval differences. |
| `pnn50` | percent | Percentage of successive interval differences above 50 ms. |

## Notes

`tdmetrics` is conformant with the Biosiglib `hrv.tdmetrics` specification. The function is modality-generic and can be used with beat or pulse interval series when the input contract is satisfied.

## Example

```python
import numpy as np

from biosigpy.hrv import tdmetrics

intervals = np.array([0.82, 0.80, 0.84, np.nan, 0.81, 0.83])
metrics = tdmetrics(intervals)

print(metrics["mhr"])
print(metrics["rmssd"])
```
