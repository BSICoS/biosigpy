# Code Style Guide

- Use Python type hints where they clarify public interfaces or reusable helpers.
- Use NumPy arrays for numeric data.
- Reuse validation helpers in `biosigpy.tools._validation` where appropriate.
- Reuse shared NaN-processing helpers in `biosigpy.tools._nan_processing` for NaN-aware filtering behavior.
- Prefer clear `snake_case` names that match canonical Biosiglib concepts where practical.
