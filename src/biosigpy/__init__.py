"""Python implementations of Biosiglib biomedical signal algorithms."""

from biosigpy._version import __version__
from biosigpy.ecg.baselineremove import (
    BaselineRemoveResult,
    BaselineRemoveWarning,
    baselineremove,
)
from biosigpy.ecg.pantompkins import pantompkins
from biosigpy.ecg.sloperange import SlopeRangeResult, sloperange
from biosigpy.hrv.fillgaps import FillGapsResult, fillgaps
from biosigpy.hrv.fdmetrics import (
    FdMetricsResult,
    FdMetricsWarning,
    SeparatedFdMetricsResult,
    fdmetrics,
)
from biosigpy.hrv.ipfm import IpfmResult, ipfm
from biosigpy.hrv.osp import OspResult, osp
from biosigpy.hrv.removefp import removefp
from biosigpy.hrv.tdmetrics import tdmetrics

__all__ = [
    "BaselineRemoveResult",
    "BaselineRemoveWarning",
    "FillGapsResult",
    "FdMetricsResult",
    "FdMetricsWarning",
    "IpfmResult",
    "OspResult",
    "SeparatedFdMetricsResult",
    "SlopeRangeResult",
    "__version__",
    "baselineremove",
    "fdmetrics",
    "fillgaps",
    "ipfm",
    "osp",
    "pantompkins",
    "removefp",
    "sloperange",
    "tdmetrics",
]
