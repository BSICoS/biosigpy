# biosigpy - Biomedical Signal Processing Library for Python

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-GPL--3.0-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/Tests-Passing-brightgreen)](tests/)
[![Version](https://img.shields.io/badge/Version-0.0.0-orange)](pyproject.toml)

Python implementation of the language-independent [Biosiglib](https://github.com/BSICoS/biosiglib) specifications for biomedical signal processing.

---

**Developed by**: [BSICoS Research Group](https://bsicos.i3a.es/)

**Status**: Active Development

## Installation

Biosigpy is currently installed from the repository source tree:

```bash
git clone https://github.com/BSICoS/biosigpy.git
cd biosigpy
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

See [Installation](docs/getting-started/installation.md) for Windows PowerShell activation, development extras, and documentation build setup.

## Biosiglib conformance

Biosigpy implements conformant Python APIs for the Biosiglib algorithms currently in scope. The root `conformance.json` pins the exact Biosiglib revision used by shared conformance tests.

See [Conformance](docs/conformance.md) for validation commands and local checkout details.

## Documentation

> **Documentation site**
>
> Visit: [https://bsicos.github.io/biosigpy/](https://bsicos.github.io/biosigpy/)

API reference pages are generated from public Python docstrings using `mkdocstrings`. Detailed installation, examples, conformance, and contribution notes live in [`docs/`](docs/).

## Support

- Report issues on [GitHub Issues](https://github.com/BSICoS/biosigpy/issues).
- Contact the development team for additional support.

## License

Biosigpy is distributed under the GNU General Public License version 3. See [LICENSE](LICENSE) for the complete license text.
