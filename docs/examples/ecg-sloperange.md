# ECG Slope-Range Respiration

Script: `examples/ecg/sloperange_example.py`

This example loads the shared Medicom MTD ECG, respiration, and R-wave timing
fixtures; derives the ECG signal with the same low-pass differentiator used by
Biosigmat; computes slope-range ECG-derived respiration; and compares the
aligned EDR series with the device respiration channel. The derivative ECG
panel overlays the masked upslope and downslope vectors plus the selected
maximum and minimum positions, so every accepted analysis window can be
inspected visually.

Run it from the repository root:

```bash
python examples/ecg/sloperange_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/ecg/sloperange_example.py --save-figure sloperange.png
```
