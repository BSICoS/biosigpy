# Installation

## Requirements

- Python 3.10–3.13
- NumPy 1.26.4 or newer
- SciPy 1.11.4 or newer

## Install from PyPI

Install the latest release and its required dependencies directly from PyPI:

```bash
python -m pip install biosigpy
```

Beta releases are awaiting validation in everyday use. To include beta releases when updating, use `python -m pip install --upgrade --pre biosigpy`; ordinary installs select stable releases by default.

## Install a downloaded release artifact

The [latest GitHub release](https://github.com/BSICoS/biosigpy/releases/latest) provides the same verified wheel and source distribution published to PyPI. Install a downloaded artifact when you need an offline or explicitly archived copy:

```bash
python -m pip install path/to/biosigpy-*.whl
```

## Verify the installation

```python
from importlib.metadata import version

import biosigpy

print(version("biosigpy"))
print(biosigpy.__version__)
```

## Get help

If installation fails, open an issue in the [Biosigpy repository](https://github.com/BSICoS/biosigpy/issues) and include your Python version, operating system, and the complete error message.
