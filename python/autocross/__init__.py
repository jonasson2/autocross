"""Simulation and correlation inference for autocorrelated time series."""

from ._testcase import testcase
from ._t3 import (
    AutocrossError,
    AutocrossWarning,
    CalibrationWarning,
    DegenerateBootstrapWarning,
    PersistenceWarning,
    T3Diagnostics,
    T3Result,
    T3Settings,
    pearson_t3,
)

__all__ = [
    "testcase", "pearson_t3", "T3Settings", "T3Result", "T3Diagnostics",
    "AutocrossError", "AutocrossWarning", "PersistenceWarning",
    "CalibrationWarning", "DegenerateBootstrapWarning",
]
