# HRV Heart-timing Reconstruction

Script: `examples/hrv/ipfm_example.py`

This example follows the Biosigmat IPFM workflow using the first 100 event
times from the shared Medicom MTD timing fixture. It reconstructs
instantaneous heart rate at 4 Hz, computes the optional TVIPFM modulating
signal, and displays both sampled outputs on their canonical time grid.

Run it from the repository root:

```bash
python examples/hrv/ipfm_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/hrv/ipfm_example.py --save-figure hrv-ipfm.png
```

The upper panel shows instantaneous heart rate in hertz. The lower panel shows
the dimensionless TVIPFM modulating signal after removing and normalizing by
the fixed 0.03 Hz mean-rate component.
