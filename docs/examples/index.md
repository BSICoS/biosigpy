# Examples

Executable examples live in the repository under `examples/`.

Run examples from the repository root after installing the `examples` optional dependency group:

```bash
python -m pip install -e ".[examples]"
```

## Available examples

- [HRV time-domain metrics](hrv-tdmetrics.md)
- [ECG Pan-Tompkins R-wave detection](ecg-pantompkins.md)

Examples that generate figures open an interactive Matplotlib window when an
interactive backend is available. In non-interactive environments, pass
`--save-figure PATH` to save the plot instead. To check the active backend:

```bash
python -c "import matplotlib; print(matplotlib.get_backend())"
```

In VS Code, open an example file, use **Run and Debug**, and select
`Biosigpy: current file with QtAgg`. The launch configuration runs the
currently open file because it uses `"program": "${file}"`.
