# HRV Missing-event Reconstruction

Script: `examples/hrv/fillgaps_example.py`

This example removes true events from a regular reference series and applies
`fillgaps` to the already false-positive-free timestamps. It compares the
reference intervals with the reconstructed intervals. Plotting remains outside
the numerical function, so normal calls to `fillgaps` have no GUI side effects.

Run it from the repository root:

```bash
python examples/hrv/fillgaps_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/hrv/fillgaps_example.py --save-figure hrv-fillgaps.png
```
