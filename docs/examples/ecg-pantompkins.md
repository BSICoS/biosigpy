# ECG Pan-Tompkins R-Wave Detection

Script: `examples/ecg/pantompkins_example.py`

This example loads the shared Medicom MTD ECG fixture, detects R waves with `biosigpy.ecg.pantompkins`, prints the canonical outputs, and plots the processing stages.

Run it from the repository root:

```bash
python examples/ecg/pantompkins_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/ecg/pantompkins_example.py --save-figure pantompkins.png
```
