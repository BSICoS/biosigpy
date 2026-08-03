# Installation

## Supported runtime

Biosigpy 0.1.0 supports Python 3.10 through 3.13, NumPy 1.26.4 or newer, and
SciPy 1.11.4 or newer. CI runs the complete conformance suite on every supported
Python version and also tests the minimum NumPy/SciPy pair on Python 3.10.

## Install a release artifact

Download the wheel or source archive and `SHA256SUMS` from the
[v0.1.0 GitHub release](https://github.com/BSICoS/biosigpy/releases/tag/v0.1.0).
After checking the recorded hash, install either artifact in a virtual
environment:

```bash
python -m pip install biosigpy-0.1.0-py3-none-any.whl
```

The source archive is installed in the same way:

```bash
python -m pip install biosigpy-0.1.0.tar.gz
```

Confirm the installed implementation version programmatically:

```python
from importlib.metadata import version

import biosigpy

assert version("biosigpy") == biosigpy.__version__ == "0.1.0"
```

## Development checkout

### Clone the repository

```bash
git clone https://github.com/BSICoS/biosigpy.git
cd biosigpy
```

### Create a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

### Install optional dependency groups

Install the examples extra when you want to run the executable examples:

```bash
python -m pip install -e ".[examples]"
```

Install the development extra when you want to run tests:

```bash
python -m pip install -e ".[dev]"
```

Install the documentation extra when you want to build or serve the docs:

```bash
python -m pip install -e ".[docs]"
```

### Build or serve the documentation

```bash
mkdocs build --strict
mkdocs serve
```
