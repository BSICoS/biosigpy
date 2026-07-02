"""Example usage of biosigpy.ecg.pantompkins."""

import numpy as np

from biosigpy.ecg import pantompkins


sampling_frequency = 256.0
duration = 10.0
time = np.arange(0.0, duration, 1.0 / sampling_frequency)

# Build a simple synthetic ECG-like signal with narrow pulses.
rng = np.random.default_rng(42)
ecg = 0.03 * rng.standard_normal(time.size)
ecg += 0.05 * np.sin(2.0 * np.pi * 0.3 * time)

reference_times = np.arange(0.8, duration, 1.0)
for reference_time in reference_times:
    ecg += np.exp(-0.5 * ((time - reference_time) / 0.015) ** 2)

outputs = pantompkins(ecg, sampling_frequency)

print("Detected R-wave times in seconds:")
print(outputs["r_wave_times"])
print("Output arrays:")
for output_name, output_value in outputs.items():
    print(f"{output_name}: shape={output_value.shape}")
