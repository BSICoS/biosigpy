# HRV False-positive Removal

Script: `examples/hrv/removefp_example.py`

This example loads the first 50 event times from the same
`medicom_mtd_r_wave_timing.csv` fixture as Biosigmat. It inserts false-positive
events after MATLAB beats 10, 20, and 30 with the same 0.05, 0.08, and 0.06 s
offsets, then applies `removefp`. It plots the event timestamps and interval
series before and after correction. The numerical function has no GUI side
effects; plotting belongs only to this example.

Run it from the repository root:

```bash
python examples/hrv/removefp_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/hrv/removefp_example.py --save-figure hrv-removefp.png
```
