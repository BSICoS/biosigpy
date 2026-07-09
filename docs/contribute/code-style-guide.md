# Code Style Guide

- Use Python type hints where they clarify public interfaces or reusable helpers.
- Use NumPy arrays for numeric data.
- Reuse validation helpers in `biosigpy.tools._validation` where appropriate.
- Reuse shared NaN-processing helpers in `biosigpy.tools._nan_processing` for NaN-aware filtering behavior.
- When porting algorithm logic, keep variable names close to Biosigmat or MATLAB names converted to Python `snake_case`.
