# AGENTS

Persistent project rules for coding agents working in Biosigpy:

1. Biosigpy is the Python implementation of the Biosiglib contracts. The pinned Biosiglib JSON specifications and shared cases define normative behavior.
2. Preserve formulas, units, defaults, filtering direction and phase, NaN behavior, edge cases, physiological meaning, and reference results. Ask the maintainer before changing scientific behavior.
3. Use idiomatic Python APIs and internals when they do not alter the contract. Canonical public IDs normally remain `snake_case`; ECG timing uses `r_wave_*`.
4. Keep generic interval methods modality-neutral unless their contract explicitly narrows them.
5. `biosiglib.lock` contains one exact lowercase Biosiglib commit. Merged code must pass every specification and shared case in that commit; partial work belongs in issues or pull requests.
6. Consume fixtures and cases from the Biosiglib checkout. Do not copy them into this repository.
7. The normal full pytest run must validate the lock, verify the checkout commit, and execute all discovered shared cases.
8. Put executable examples under `examples/` and describe what they teach a user, not their implementation ancestry.
9. Use English for code, comments, filenames, and technical documentation. Use the repository-local `.venv`; do not commit environments, caches, editor settings, or downloaded data.
10. Avoid generic resource APIs and unnecessary cross-language infrastructure.
