# Installation

Biosigpy is currently installed from the repository source tree.

## Clone the repository

```bash
git clone https://github.com/BSICoS/biosigpy.git
cd biosigpy
```

## Create a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

## Install optional dependency groups

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

## Build or serve the documentation

```bash
mkdocs build --strict
mkdocs serve
```
