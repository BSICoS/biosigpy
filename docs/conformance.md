# Conformance

Biosigpy uses `conformance.json` to pin the exact Biosiglib revision used for shared conformance testing. The manifest declares conformance with every specification in that revision; it is not a support inventory or roadmap.

Shared conformance cases define scientific and computational behavior across implementations. Python-specific tests additionally cover Python API behavior, packaging, exceptions, and local implementation details.

The test suite discovers every specification and case file from the pinned Biosiglib checkout. Case definitions determine whether the implementation is expected to return outputs or raise an implementation-specific error. A full default pytest run fails during collection if a specification has no shared cases or any discovered case is not collected, so a future Biosiglib commit cannot silently extend the contract without exercising its new coverage.

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
