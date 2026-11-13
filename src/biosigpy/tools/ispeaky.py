"""Classify spectra using both peakedness thresholds."""

import numpy as np
from numpy.typing import ArrayLike


def ispeaky(pkl: ArrayLike, akl: ArrayLike, pkl_threshold: float,
            akl_threshold: float) -> bool | np.ndarray:
    """Return whether both peakedness measures meet their inclusive thresholds.

    PKL and AKL must be nonempty real numeric scalars or arrays of the same
    shape. No broadcasting is performed. Each threshold must be a finite
    scalar between 0 and 100. A NaN measure produces False at that position.
    Array inputs retain their shape; scalar inputs return a bool.
    """
    # Parse and validate inputs
    pkl_values = np.asarray(pkl)
    akl_values = np.asarray(akl)
    for name, values in (("pkl", pkl_values), ("akl", akl_values)):
        if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
            raise TypeError(f"{name} must be real numeric data")
        if values.size == 0:
            raise ValueError(f"{name} must be nonempty")
    for name, threshold in (("pkl_threshold", pkl_threshold), ("akl_threshold", akl_threshold)):
        value = np.asarray(threshold)
        if not np.issubdtype(value.dtype, np.number) or np.iscomplexobj(value):
            raise TypeError(f"{name} must be real numeric data")
        if value.ndim != 0 or not np.isfinite(value) or not 0 <= value <= 100:
            raise ValueError(f"{name} must be a finite scalar between 0 and 100")

    # Validate that pkl and akl have the same dimensions
    if pkl_values.shape != akl_values.shape:
        raise ValueError("pkl and akl must have the same shape")

    # Determine if spectra are peaky based on both criteria
    is_peaky = (pkl_values >= pkl_threshold) & (akl_values >= akl_threshold)
    return bool(is_peaky) if is_peaky.ndim == 0 else is_peaky
