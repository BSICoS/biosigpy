# Biosigpy

**Python implementation of the Biosiglib specifications for biomedical signal processing.**

Biosigpy provides Python implementations of the language-independent algorithms and behavior defined by [Biosiglib](https://github.com/BSICoS/biosiglib). It is part of the BSICoS biomedical signal-processing ecosystem together with [Biosigmat](https://github.com/BSICoS/biosigmat), the MATLAB implementation.

## Status

Biosigpy is in its initial development stage and is not yet available as a stable package.

## Installation for development

Clone the repository and run examples or tests from the source tree:

```bash
git clone https://github.com/BSICoS/biosigpy.git
cd biosigpy
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy pytest
```

On Windows PowerShell, use `.venv\Scripts\Activate.ps1` instead of `source .venv/bin/activate`.

Until packaging is finalized, run scripts with `PYTHONPATH=src` from the repository root.

## Basic usage

```python
import numpy as np

from biosigpy.hrv import tdmetrics

intervals = np.array([0.82, 0.80, 0.84, np.nan, 0.81, 0.83])
metrics = tdmetrics(intervals)

print(metrics["mhr"])
print(metrics["rmssd"])
```

More examples are available under [`examples/`](examples/):

```bash
PYTHONPATH=src python examples/hrv/tdmetrics_example.py
PYTHONPATH=src python examples/ecg/pantompkins_example.py
```

## Documentation

The Biosigpy documentation site is being developed. The planned location is [https://bsicos.github.io/biosigpy/](https://bsicos.github.io/biosigpy/).

Language-independent algorithm behavior is documented centrally in [Biosiglib](https://github.com/BSICoS/biosiglib). Biosigpy documentation should focus on Python installation, API usage, executable examples, and Python-specific implementation notes.

## Conformance

Biosigpy includes `conformance.json`, a machine-readable manifest that pins the exact Biosiglib revision used for shared conformance tests.

Shared conformance cases define the scientific and computational behavior that must remain aligned across implementations. Python-specific tests may additionally cover Python API behavior, data types, exceptions, packaging, and internal implementation details.

## Contributing

Changes that alter scientific or computational behavior must be reflected in Biosiglib and reviewed explicitly before implementation. Language-specific choices are fine when they do not change the normative behavior.

## License

Biosigpy is distributed under the GNU General Public License version 3. See [LICENSE](LICENSE) for the complete license text.
