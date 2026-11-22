"""Public interface to the uncorrected PearsonT3 interval."""

from dataclasses import dataclass
from numbers import Real
import operator
from threading import Lock
import warnings

import numpy as np
import randompack

from ._core import run_t3

_ENGINE_LOCK = Lock()
_C_INT_MAX = np.iinfo(np.intc).max


class AutocrossError(RuntimeError):
    """PearsonT3 could not compute a usable result."""


class AutocrossWarning(UserWarning):
    """Base class for PearsonT3 warnings."""


class PersistenceWarning(AutocrossWarning):
    """Persistence estimation used a boundary value or a bounded refit."""


class CalibrationWarning(AutocrossWarning):
    """The calibration grid does not bracket the requested confidence."""


class DegenerateBootstrapWarning(AutocrossWarning):
    """A bootstrap resample was constant or had no inner-bootstrap spread."""


def _count(value, name, minimum):
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if not minimum <= value <= _C_INT_MAX:
        raise ValueError(f"{name} must be in [{minimum}, {_C_INT_MAX}]")
    return value


@dataclass(frozen=True, kw_only=True)
class T3Settings:
    """Reusable PearsonT3 settings.

    Parameters
    ----------
    b1, b2 : int, optional
        First-level resamples and second-level resamples per first-level
        resample, respectively. Defaults are 2000 and 1000.
    n_lambda : int, optional
        Number of points on the calibration grid. Default is 500.
    block_length : int or None, optional
        Block length in observations. None selects the classical rule
        based on the estimated persistence times.
    detrend : {"mean", "linear"}, optional
        Remove the mean (default) or a fitted linear trend from each series.
    """

    b1: int = 2000
    b2: int = 1000
    n_lambda: int = 500
    block_length: int | None = None
    detrend: str = "mean"

    def __post_init__(self):
        for name, minimum in (("b1", 2), ("b2", 2), ("n_lambda", 1)):
            value = _count(getattr(self, name), name, minimum)
            object.__setattr__(self, name, value)
        if self.block_length is not None:
            object.__setattr__(self, "block_length", _count(self.block_length, "block_length", 1))
        if not isinstance(self.detrend, str):
            raise TypeError("detrend must be a string")
        if self.detrend not in ("mean", "linear"):
            raise ValueError("detrend must be 'mean' or 'linear'")


@dataclass(frozen=True, eq=False)
class T3Diagnostics:
    """Intermediate quantities from the classical calibration.

    se_z is the first-level bootstrap SD on the Fisher scale.
    bootstrap_z and inner_se_z hold the first-level Fisher-transformed
    correlations and their second-level SE estimates, in resample order.
    calibration_curve gives the estimated coverage at each lambda_grid point.
    selected_lambda is the lambda used to construct the final interval.
    persistence_refit flags a bounded lag-one refit for x and y, respectively.
    """

    se_z: float
    selected_lambda: float
    bootstrap_z: np.ndarray
    inner_se_z: np.ndarray
    lambda_grid: np.ndarray
    calibration_curve: np.ndarray
    persistence_refit: tuple[bool, bool] = (False, False)


@dataclass(frozen=True)
class T3Result:
    """Pearson correlation and confidence interval for the detrended series.

    ci is a (lower, upper) tuple in correlation units. taux and tauy are
    estimated persistence times in the units of t. block_length is the
    actual block length, in observations. diagnostics is None unless requested.
    """

    r: float
    ci: tuple[float, float]
    taux: float
    tauy: float
    block_length: int
    confidence: float
    settings: T3Settings
    diagnostics: T3Diagnostics | None = None


def _vector(value, name):
    values = np.asarray(value)
    if values.dtype.kind not in "iuf":
        raise TypeError(f"{name} must be real numeric data")
    if values.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    with np.errstate(over="ignore"):
        values = np.ascontiguousarray(values, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError(f"{name} must contain only finite values")
    return values


def pearson_t3(t, x, y, *, confidence=0.95, settings=None, rng=None, diagnostics=False):
    """Compute the uncorrected PearsonT3 correlation confidence interval.

    Parameters
    ----------
    t, x, y : array_like
        Observation times and paired observations. Each input must be a
        finite, real, one-dimensional array of the same length, at least
        10. Times must be strictly increasing.
    confidence : float, optional
        Two-sided confidence level in (0, 1). Default is 0.95.
    settings : T3Settings, optional
        Bootstrap and preprocessing settings. None uses the defaults.
    rng : randompack.Rng, optional
        Random number generator for the bootstrap. A supplied generator
        advances in place; otherwise a temporary generator is used.
    diagnostics : bool, optional
        Return intermediate bootstrap and calibration quantities in
        result.diagnostics. Default is False.

    Returns
    -------
    T3Result
        Correlation of the detrended series, confidence interval,
        estimated persistence times, and block length.

    Raises
    ------
    ValueError
        The data or settings are invalid, or a series is constant after
        detrending.
    AutocrossError
        The numerical calculation cannot produce a usable result.

    Notes
    -----
    Uses the paired moving-block bootstrap, Fisher transformation, and
    classical PearsonT3 calibration. Warns if a persistence bound or refit is used,
    available bootstrap coverage does not bracket the requested level, or
    some resamples are degenerate but a usable interval is returned.

    Examples
    --------
    >>> from autocross import T3Settings, pearson_t3, testcase
    >>> from randompack import Rng
    >>> rng = Rng()
    >>> rng.seed(123)
    >>> t, x, y = testcase(80, model="ou", spacing=("gamma", 16), rho=0.3, rng=rng)
    >>> result = pearson_t3(t, x, y, settings=T3Settings(b1=400, b2=200), rng=rng)
    >>> lower, upper = result.ci
    >>> lower <= result.r <= upper
    True
    """
    return _pearson_t3_levels(t, x, y, levels=(confidence,), settings=settings, rng=rng,
                              diagnostics=diagnostics)[0]


def _pearson_t3_levels(t, x, y, *, levels, settings=None, rng=None, diagnostics=False):
    """Evaluate multiple confidence levels on a single bootstrap curve."""
    if settings is None:
        settings = T3Settings()
    elif not isinstance(settings, T3Settings):
        raise TypeError("settings must be a T3Settings instance")
    levels = tuple(levels)
    if not levels:
        raise ValueError("at least one confidence level is required")
    for confidence in levels:
        if not isinstance(confidence, Real) or isinstance(confidence, (bool, np.bool_)):
            raise TypeError("confidence must be a real number")
        if not 0 < float(confidence) < 1:
            raise ValueError("confidence must be in (0, 1)")
    levels = np.array(levels, dtype=np.float64)
    if not isinstance(diagnostics, (bool, np.bool_)):
        raise TypeError("diagnostics must be a boolean")
    if rng is None:
        rng = randompack.Rng()
    elif not isinstance(rng, randompack.Rng):
        raise TypeError("rng must be a randompack.Rng")

    t, x, y = (_vector(v, name) for v, name in zip((t, x, y), ("t", "x", "y")))
    n = t.size
    if x.size != n or y.size != n:
        raise ValueError("t, x, and y must have the same length")
    if not 10 <= n <= _C_INT_MAX // 2:
        raise ValueError(f"number of observations must be in [10, {_C_INT_MAX // 2}]")
    if not (t[1:] > t[:-1]).all() or not np.isfinite(float(t[-1]) - float(t[0])):
        raise ValueError("t must be strictly increasing with a finite time span")
    if (x == x[0]).all() or (y == y[0]).all():
        raise ValueError("x and y must not be constant")
    block = settings.block_length
    if block is not None and block >= n:
        raise ValueError("block_length must be smaller than the number of observations")

    with _ENGINE_LOCK:
        status, output = run_t3(
            t, x, y, levels, settings.b1, settings.b2, settings.n_lambda,
            block or 0, settings.detrend == "linear", rng,
        )
    if status == 1:
        raise MemoryError("could not allocate PearsonT3 workspace")
    if status == 2:
        failed = " and ".join(name for name, state in zip(("x", "y"), output) if state < 0)
        raise AutocrossError(f"could not estimate persistence for {failed}")
    if status == 5:
        raise ValueError("x and y must not be constant after detrending")
    if status:
        raise AutocrossError("PearsonT3 could not compute a finite interval with usable bootstrap spread")

    r, ci, (taux, tauy), actual, se_z, selected, cover, z1, se2, grid, states, constant, zero = output
    refit = tuple(bool(state & 2) for state in states)
    for bit, outcome in ((1, "a boundary persistence estimate"), (2, "a bounded persistence refit")):
        names = " and ".join(name for name, state in zip(("x", "y"), states) if state & bit)
        if names:
            warnings.warn(f"PearsonT3 used {outcome} for {names}", PersistenceWarning, stacklevel=3)
    low, high = cover.min(), cover.max()
    if constant or zero:
        warnings.warn(
            f"PearsonT3 encountered {constant} constant resamples and "
            f"{zero} first-level samples with zero inner-bootstrap SE",
            DegenerateBootstrapWarning, stacklevel=3,
        )
    if diagnostics:
        for array in (z1, se2, grid, cover):
            array.setflags(write=False)
    results = []
    for confidence, ends, choice in zip(levels, ci, selected):
        if not low <= confidence <= high:
            warnings.warn(
                f"requested confidence {confidence:g} is outside the calibration coverage range "
                f"[{low:g}, {high:g}]; returning the closest available interval",
                CalibrationWarning, stacklevel=3,
            )
        info = T3Diagnostics(se_z, float(choice), z1, se2, grid, cover, refit) if diagnostics else None
        results.append(T3Result(r, tuple(map(float, ends)), taux, tauy, actual,
                                float(confidence), settings, info))
    return results
