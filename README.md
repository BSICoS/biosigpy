# Biosigpy

**Python implementation of the Biosiglib specifications for biomedical signal processing.**

Biosigpy is the Python implementation of the language-independent algorithms and behavior defined by [Biosiglib](../biosiglib).

It is part of an ecosystem that currently includes:

* [Biosiglib](../biosiglib): specifications, scientific provenance, shared fixtures, conformance cases, and common usage scenarios.
* [Biosigmat](../biosigmat): MATLAB implementation.
* **Biosigpy**: Python implementation.

## Purpose

Biosigpy provides reproducible biomedical signal-processing tools for multiple signal modalities, including:

* ECG.
* PPG.
* Respiration.
* Heart rate variability.
* Related biomedical signal-processing utilities.

The project is intended to reproduce the scientific and computational behavior defined by Biosiglib while providing an idiomatic Python interface.

## Independent implementation

Biosigpy is not a line-by-line translation or clone of Biosigmat.

Both libraries implement the same Biosiglib specifications, but they may differ in:

* Internal architecture.
* Function and module organisation.
* Python and MATLAB naming conventions.
* Data structures.
* Error types.
* Plotting libraries.
* Performance optimisations.

They must remain consistent in the aspects defined normatively by Biosiglib, including:

* Scientific meaning.
* Input and output units.
* Default parameters.
* Mathematical and computational behavior.
* Missing-value handling.
* Edge-case behavior.
* Numerical results within the established tolerances.

## Development status

Biosigpy is in its initial development stage and is not yet available as a stable package.

The first implemented algorithms will be:

* `hrv.tdmetrics`
* `ecg.pantompkins`

These pilots will be used to validate the Biosiglib architecture with:

* A simple numerical function.
* A more complex signal-processing algorithm.
* Shared fixtures.
* Cross-language conformance cases.
* Automated comparison with Biosigmat.

The rest of the public Biosigmat API will be incorporated progressively after the pilot architecture has been validated.

## Specifications

Biosigpy implements machine-readable specifications maintained in Biosiglib.

Each specification defines:

* Inputs and outputs.
* Units.
* Parameters and default values.
* Mathematical or computational behavior.
* Missing-value behavior.
* Edge cases.
* Absolute numerical tolerances.
* Scientific references.
* Associated conformance cases.

Biosigpy does not redefine this information independently.

## Conformance testing

Biosigpy uses its native Python testing tools to run the shared Biosiglib conformance cases.

The shared resources define:

* Input fixtures.
* Parameters.
* Expected outputs.
* Absolute comparison tolerances.
* Expected `NaN` behavior.

Python-specific unit tests may additionally verify:

* Public API behavior.
* Python data types.
* Exceptions.
* Package integration.
* Internal implementation details.

These implementation-specific tests complement the shared conformance suite but do not replace it.

## Shared examples and workflows

Biosigpy will provide executable Python versions of the common examples and workflows described by Biosiglib.

Corresponding examples in Biosigpy and Biosigmat should follow the same scientific workflow, for example:

1. Load the same signal or an equivalent fixture.
2. Apply the same preprocessing stages.
3. Use equivalent algorithm parameters.
4. Detect or calculate the same physiological events or measurements.
5. Present comparable results and plots.

The Python code and visualisation libraries may differ from MATLAB, but the example should teach the same scientific process and lead to equivalent interpretations.

## Conformance manifest

Biosigpy will include a machine-readable conformance manifest identifying:

* The Biosigpy version.
* The Biosiglib release and commit used for validation.
* The specifications currently implemented.
* Their conformance status.

Continuous integration will verify these declarations against the referenced Biosiglib version.

## Documentation

Biosigpy documentation will include:

* Installation instructions.
* Python API reference.
* Executable examples.
* Development guidance.
* Implementation-specific behavior.

Scientific definitions and language-independent behavior will be documented centrally in Biosiglib to avoid maintaining conflicting copies.

## Installation

Biosigpy has not yet published its first installable release.

Installation and development-environment instructions will be added when the initial package structure and supported Python versions have been established.

## Contributing

Contributions should preserve the separation between:

* Language-independent behavior defined by Biosiglib.
* Python-specific implementation decisions defined in Biosigpy.

Changes that alter scientific or computational behavior must first be reflected in the corresponding Biosiglib specification.

## License

Biosigpy is distributed under the GNU General Public License version 3.

See `LICENSE` for the complete license text.
