# HRV Missing-event Reconstruction

Script: `examples/hrv/fillgaps_example.py`

This example loads the first 100 event times from the same
`medicom_mtd_r_wave_timing.csv` fixture as Biosigmat. It removes the exact
detections selected by MATLAB with `rng(40)`, applies the explicit `removefp`
then `fillgaps` sequence, and compares the reference intervals with the
reconstructed intervals. Normal calls to `fillgaps` have no GUI side effects.

Run it from the repository root:

```bash
python examples/hrv/fillgaps_example.py
```

To inspect the reconstruction interactively, run:

```bash
python examples/hrv/fillgaps_example.py --debug
```

The upper panel shows the current intervals, detected gaps, and adaptive
detection threshold. The lower panel shows each attempted reconstruction in
green or red together with its validation limits. Processing waits for a key
press or mouse click after every attempt and closes the debug figure when it
finishes. This mode requires Matplotlib with an interactive backend.

To save the figure instead of opening an interactive window:

```bash
python examples/hrv/fillgaps_example.py --save-figure hrv-fillgaps.png
```
