# AGENTS.md

This file provides persistent context for AI coding agents working in Biosigpy. Treat the decisions below as project policy unless the maintainers explicitly change them.

- Biosigpy is the Python implementation of the language-independent Biosiglib specifications.
- Biosiglib is the source of truth for normative public algorithm behavior, units, defaults, missing-value handling, edge cases, scientific provenance, shared fixtures, and conformance cases.
- Biosigpy is not a line-by-line translation of Biosigmat. It may use idiomatic Python APIs, module organization, data structures, exceptions, and plotting libraries while preserving normative behavior.
- Biosigmat is the mature starting implementation and may provide reference outputs for complex algorithms. It is not automatically correct when a disagreement is detected.
- If Biosigpy and Biosigmat disagree in a scientifically meaningful way, document the disagreement and ask the maintainer before changing Biosigpy, Biosigmat, or the Biosiglib specification.
- Do not change scientific or computational behavior without explicit maintainer review. This includes filtering direction and phase behavior, NaN handling, default filters, default parameters, units, physiological interpretation, and reference-result provenance.
- A change from zero-phase NaN-aware filtering to causal filtering, a change from `nanfiltfilt`-style behavior to `lfilter`-style behavior, or a replacement of the `lpdfilter`-derived default with an ad-hoc default is a scientific/computational decision. Ask the maintainer before making that change.
- Do not ask about purely idiomatic differences unless they affect scientific behavior. Examples: zero-based versus one-based internal indexing, exception class names, plotting library choices, or local variable names normally do not require maintainer escalation.
- Canonical Biosiglib IDs use snake_case. Python APIs should generally use those names directly.
- ECG R-wave timing identifiers must use `r_wave_*`, not `r_peak_*`.
- Generic timing or interval algorithms, including `hrv.tdmetrics`, must remain modality-generic unless the contract explicitly narrows them.
- Examples remain under `examples/`; corresponding examples across implementations should preserve the same scientific workflow where practical.
- `conformance.json` pins an exact Biosiglib commit.
- `conformant` may only be declared after all applicable shared cases pass.
- Shared fixtures and cases must be consumed from a Biosiglib checkout and not copied back into Biosigpy.
- Implementation manifest validation must use Biosiglib's `tools/validate_specs.py --manifest PATH` command instead of duplicating JSON Schema validation code.
- Code, comments, filenames, and technical documentation must be in English.
- Avoid generic resource APIs and unnecessary cross-language infrastructure.

## Local change workflow

- Work in a local checkout for multi-file changes.
- Make related changes locally, run validation locally where possible, and push one coherent commit or a small coherent commit series.
- Use a repository-local `.venv` for Python tooling. Invoke the `.venv` Python executable explicitly in scripts and documentation when practical.
- Do not commit `.venv`, generated Python caches, local editor settings, or downloaded external data.
