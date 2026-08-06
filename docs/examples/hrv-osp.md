# HRV Respiration-related Decomposition

Script: `examples/hrv/osp_example.py`

This example follows the Biosigmat workflow using the local Medicom MTD
fixtures. It reconstructs the TVIPFM modulating signal, aligns the detrended
respiration signal to the same 4 Hz grid, estimates its spectrum, and applies
OSP to separate respiration-related modulation from the orthogonal residual.
The four panels reproduce the Biosigmat visualization: aligned respiration,
the full modulation signal, the unrelated residual, and the
respiration-related component.

Run it from the repository root:

```bash
python examples/hrv/osp_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/hrv/osp_example.py --save-figure hrv-osp.png
```
