# Conformance

Biosigpy uses `conformance.json` to pin the exact Biosiglib revision used for shared conformance testing.

Shared conformance cases define scientific and computational behavior across implementations. Python-specific tests additionally cover Python API behavior, packaging, exceptions, and local implementation details.

Scientific or computational behavior changes should first be reflected in Biosiglib before implementation-specific changes are made.

## Run validation on Windows PowerShell

```powershell
$env:PYTHONPATH='src'; .venv\Scripts\python.exe -m pytest
.venv\Scripts\python.exe ..\biosiglib\tools\validate_specs.py --manifest conformance.json
```

## Run validation on Unix shells

```bash
PYTHONPATH=src python -m pytest
python ../biosiglib/tools/validate_specs.py --manifest conformance.json
```
