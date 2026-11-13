---
title: Citation
---

# Cite Biosigpy

Research using Biosigpy should cite both the exact released Biosigpy version
that was executed and the original publications for every scientific method
used.

## Identify the version

Obtain the installed version in Python:

```python
import biosigpy

print(biosigpy.__version__)
```

Use the [Biosigpy Zenodo collection](https://doi.org/10.5281/zenodo.22112727)
to select the version that was executed and cite its version-specific DOI.

GitHub provides citation formats generated from the project's structured
metadata.

## Cite the methods

Find the original publications for each method in the
[Biosiglib method catalog](https://bsicos.github.io/biosiglib/methods/).
Biosiglib is the language-independent contract and reference catalog, not the
executable software citation.

## Unreleased code

If the analysis used an unreleased commit, cite the
[Biosigpy repository](https://github.com/BSICoS/biosigpy) together with the full
commit SHA instead of a released version.
