"""Electrocardiogram signal-processing algorithms."""

from biosigpy.ecg.pantompkins import pantompkins
from biosigpy.ecg.sloperange import SlopeRangeResult, sloperange

__all__ = ["SlopeRangeResult", "pantompkins", "sloperange"]
