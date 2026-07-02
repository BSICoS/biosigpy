"""Example usage of biosigpy.hrv.tdmetrics."""

import numpy as np

from biosigpy.hrv import tdmetrics


# Cleaned beat-to-beat or pulse-to-pulse intervals in seconds.
# NaN values are allowed as missing or invalid interval markers.
intervals = np.array([0.82, 0.80, 0.84, np.nan, 0.81, 0.83], dtype=float)

metrics = tdmetrics(intervals)

print("Time-domain variability metrics")
for metric_name, metric_value in metrics.items():
    print(f"{metric_name}: {metric_value:.3f}")
