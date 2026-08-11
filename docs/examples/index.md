# Examples

Executable examples live in the repository under `examples/`.

Run examples from the repository root after installing the `examples` optional dependency group:

```bash
python -m pip install -e ".[examples]"
```

## Available examples

- [HRV false-positive removal](hrv-removefp.md)
- [HRV missing-event reconstruction](hrv-fillgaps.md)
- [HRV frequency-domain metrics](hrv-fdmetrics.md)
- [HRV heart-timing reconstruction](hrv-ipfm.md)
- [HRV respiration-related decomposition](hrv-osp.md)
- [HRV time-domain metrics](hrv-tdmetrics.md)
- [ECG baseline-wander removal](ecg-baselineremove.md)
- [ECG Pan-Tompkins R-wave detection](ecg-pantompkins.md)
- [ECG slope-range respiration](ecg-sloperange.md)

Examples that generate figures open an interactive Matplotlib window when an
interactive backend is available. In non-interactive environments, pass
`--save-figure PATH` to save the plot instead. To check the active backend:

```bash
python -c "import matplotlib; print(matplotlib.get_backend())"
```

In VS Code, open an example file, use **Run and Debug**, and select
`Biosigpy: current file with QtAgg`. The launch configuration runs the
currently open file because it uses `"program": "${file}"`.
