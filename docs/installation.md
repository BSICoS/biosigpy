# Installation

## Requirements

- Python 3.10–3.13
- NumPy 1.26.4 or newer
- SciPy 1.11.4 or newer

## Install the latest release

Download a wheel from the [latest GitHub release](https://github.com/BSICoS/biosigpy/releases/latest), then install the downloaded file in your environment:

```bash
python -m pip install path/to/biosigpy-*.whl
```

The source archive available on the same release page can be installed in the same way.

## Verify the installation

```python
from importlib.metadata import version

import biosigpy

assert version("biosigpy") == biosigpy.__version__
print(biosigpy.__version__)
```

## Get help

If installation fails, open an issue in the [Biosigpy repository](https://github.com/BSICoS/biosigpy/issues) and include your Python version, operating system, and the complete error message.
