"""Electrocardiogram signal-processing algorithms."""

from biosigpy.ecg.baselineremove import (
    BaselineRemoveResult,
    BaselineRemoveWarning,
    baselineremove,
)
from biosigpy.ecg.pantompkins import pantompkins
from biosigpy.ecg.sloperange import SlopeRangeResult, sloperange

__all__ = [
    "BaselineRemoveResult",
    "BaselineRemoveWarning",
    "SlopeRangeResult",
    "baselineremove",
    "pantompkins",
    "sloperange",
]
