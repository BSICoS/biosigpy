# ECG Baseline-Wander Removal

Script: `examples/ecg/baselineremove_example.py`

This is the Python port of the Biosigmat `baselineremoveExample`. It loads the
same Medicom MTD ECG and R-wave timing fixtures, uses the same 256 Hz sampling
frequency, places fiducials 150 ms before each R wave, and applies the default
local averaging window. The two aligned panels reproduce the original example:
the ECG with its estimated baseline, R waves, and fiducial points, followed by
the corrected ECG with the same markers.

Run it from the repository root:

```bash
python examples/ecg/baselineremove_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/ecg/baselineremove_example.py --save-figure baselineremove.png
```
