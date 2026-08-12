# Conformance

Biosigpy uses the one-line `biosiglib.lock` file to pin the exact Biosiglib revision used for shared conformance testing. Merged code must conform to every specification in that revision; partial support and roadmaps belong in issues and pull requests.

Shared conformance cases define scientific and computational behavior across implementations. Python-specific tests additionally cover Python API behavior, packaging, exceptions, and local implementation details.

The test suite discovers every specification and case file from the pinned Biosiglib checkout. Case definitions determine whether the implementation is expected to return outputs or raise an implementation-specific error. A full default pytest run fails during collection if a specification has no shared cases or any discovered case is not collected, so a future Biosiglib commit cannot silently extend the contract without exercising its new coverage.

Scientific or computational behavior changes should first be reflected in Biosiglib before implementation-specific changes are made.

## Run validation on Windows PowerShell

```powershell
$env:PYTHONPATH='src'; .venv\Scripts\python.exe -m pytest
```

## Run validation on Unix shells

```bash
PYTHONPATH=src python -m pytest
```

This single command validates the lock, verifies the Biosiglib checkout commit, and executes every shared case.
