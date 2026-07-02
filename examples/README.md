# Examples

Executable examples for the current Biosigpy pilot algorithms.

The examples use local copies of the Medicom MTD ECG fixtures under
`examples/fixtures/ecg/`, matching the corresponding Biosigmat examples.
Install the example dependencies with an editable install:

```bash
python -m pip install -e ".[examples]"
```

Run examples from the repository root:

```bash
python examples/hrv/tdmetrics_example.py
python examples/ecg/pantompkins_example.py
```

On Windows PowerShell:

```powershell
python examples\hrv\tdmetrics_example.py
python examples\ecg\pantompkins_example.py
```

Examples that generate figures open an interactive Matplotlib window when an
interactive backend is available. In non-interactive or headless environments,
pass `--save-figure PATH` to save the plot instead:

```bash
MPLBACKEND=QtAgg python examples/ecg/pantompkins_example.py
```

```bash
MPLBACKEND=Agg python examples/ecg/pantompkins_example.py --save-figure /tmp/pantompkins_example.png
```

On Windows PowerShell:

```powershell
$env:MPLBACKEND = "Agg"
python examples\ecg\pantompkins_example.py --save-figure "$env:TEMP\pantompkins_example.png"
Remove-Item Env:\MPLBACKEND
```

For VS Code, install the example dependencies, open any example file, use
**Run and Debug**, and select `Biosigpy: current file with QtAgg`. The launch
configuration runs the currently open file because it uses `"program": "${file}"`.

To check the active Matplotlib backend:

```bash
python -c "import matplotlib; print(matplotlib.get_backend())"
```

For local validation:

```bash
python -m pip install -e ".[dev]"
python examples/hrv/tdmetrics_example.py
MPLBACKEND=Agg python examples/ecg/pantompkins_example.py --save-figure /tmp/pantompkins_example.png
python -m pytest
```
