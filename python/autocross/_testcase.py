"""Stationary test cases for correlation inference."""

import math
import operator
from numbers import Real

import numpy as np
import randompack
import varmapack


def testcase(n=None, *, model, spacing=None, t=None, rho=0.0, marginal="normal",
             tau=None, drift=None, theta=None, lag1=None, rng=None):
    """Simulate a pair of stationary time series on a common observation grid.

    Supply either ``t`` for explicit observation times, or both ``n`` and
    ``spacing`` to generate a grid.

    Parameters
    ----------
    n : int, optional
        Length of each returned vector for a generated grid.
    model : {"ou", "ar1", "slepian"}
        Latent continuous-time Gaussian process: Ornstein--Uhlenbeck
        ("ou") or a unit-window Slepian moving-average process ("slepian").
        "ar1" is a synonym for "ou".
    spacing : tuple, optional
        Grid specification ``(name, parameter)``. Generated grids start at 1.
        See Spacing rules below.
    t : array_like, optional
        Nonempty one-dimensional array of finite, strictly increasing real
        observation times. A range, list, or NumPy array is accepted.
    rho : float, default 0
        Population contemporaneous Pearson correlation between the
        returned series.
    marginal : str or tuple, default "normal"
        Marginal distribution of each returned series. All choices have
        population mean zero and variance one. See Marginal distributions
        below.
    tau : array_like, shape (2,), optional
        Finite, positive persistence times for the two components of an OU
        process with diagonal drift, in the units of ``t``. The autocorrelation
        of component ``i`` at elapsed time ``d`` is ``exp(-d / tau[i])``.
        Defaults to (1, 1) when no drift matrix is supplied.
    drift : array_like, shape (2, 2), optional
        Real OU drift matrix K, an alternative to ``tau``. The transition
        over an interval d is the matrix exponential of -d*K. The
        eigenvalues of K must have positive real parts.
    theta : float or array_like, optional
        For "slepian", supply exactly one of ``theta`` and ``lag1``.
        ``theta`` is the moving-average coefficient matrix at unit spacing.
        A scalar gives a multiple of the identity, a pair gives the
        diagonal entries, and a (2, 2) matrix allows coupled dynamics.
    lag1 : float or array_like, optional
        Autocorrelation of each latent component at a time lag of one
        unit in a diagonal Slepian model. A scalar applies to both
        components; a pair specifies them separately. Each value a
        must lie in [-0.5, 0.5] and selects the coefficient b satisfying
        ``a = b / (1 + b**2)`` and ``abs(b) <= 1``.
    rng : randompack.Rng, optional
        Random number generator for the grid and process. A supplied
        generator advances in place; otherwise use a temporary generator.

    Returns
    -------
    t, x, y : ndarray
        Observation times and the two series sampled at those times.
        Each array is one-dimensional and has dtype float64.

    Raises
    ------
    ValueError
        The grid or model parameters are invalid, or the requested
        correlation is infeasible for the selected model and marginals.

    Notes
    -----
    **Spacing rules**

    - ``("gamma", shape)``: draw Gamma gaps with a finite positive shape
      and rescale them to have realized mean one. Requires n >= 2.
    - ``("geometric", miss_prob)``: retain the first point of a unit grid,
      then omit each subsequent point independently with probability
      ``miss_prob`` in [0, 1) until n observations have been retained.
    - ``("single-gap", n_missing)``: choose an interior gap uniformly.
      Skip a positive integer number of points there and use unit
      spacing elsewhere. Give an integer for a fixed count or an
      inclusive pair of integers to draw the count independently and
      uniformly. Requires n >= 2.

    **Marginal distributions**

    - ``"normal"``: return the latent Gaussian components.
    - ``"lognormal"``: standardized lognormal with log-scale SD 1,
      corresponding to population skewness about 6.18.
    - ``("lognormal", {"skewness": g})``: standardized lognormal with
      population skewness g for both series. At g=0, return the
      Gaussian components.

    Standardization uses population moments, so standardized lognormal
    observations can be negative. The process starts in stationarity and
    continues across gaps. With diagonal dynamics, rho=0 gives independent
    series. With coupled dynamics, the series may still be correlated at
    nonzero lags.

    Examples
    --------
    >>> from autocross import testcase
    >>> t, x, y = testcase(50, model="ar1", tau=(1, 2), rho=0.3,
    ...                    spacing=("gamma", 16), marginal="lognormal")
    >>> t, x, y = testcase(model="ou", t=range(2, 22, 2))
    >>> t, x, y = testcase(50, model="slepian", lag1=0.4,
    ...                    spacing=("single-gap", (3, 8)),
    ...                    marginal=("lognormal", {"skewness": 2}))
    """
    latent_rho, sigma2, cv = _latent_parameters(rho, marginal)
    process = _make_process(model, tau, drift, theta, lag1, latent_rho)
    if rng is None:
        rng = randompack.Rng()
    elif not isinstance(rng, randompack.Rng):
        raise TypeError("rng must be a randompack.Rng")

    times = _observation_times(n, spacing, t, rng)
    try:
        values = process.sim(times, rng=rng)[0]
    except varmapack.VarmapackError as exc:
        raise ValueError(f"could not simulate at the specified times: {exc}") from exc
    if sigma2:
        with np.errstate(over="ignore"):
            values = np.expm1(math.sqrt(sigma2) * values - sigma2 / 2) / cv
    if not np.isfinite(values).all():
        raise ValueError("simulated observations are not representable as float64")
    return times, values[:, 0].copy(), values[:, 1].copy()


def _finite_real(value, name):
    if (
        not isinstance(value, Real)
        or isinstance(value, (bool, np.bool_))
        or not math.isfinite(value)
    ):
        raise ValueError(f"{name} must be a finite real number")
    return float(value)


def _real_array(value, name):
    array = np.asarray(value)
    if array.dtype.kind not in "iuf" or not np.isfinite(array).all():
        raise ValueError(f"{name} must contain finite real numbers")
    with np.errstate(over="ignore"):
        array = array.astype(np.float64, copy=False)
    if not np.isfinite(array).all():
        raise ValueError(f"{name} must be representable as float64")
    return array


def _positive_int(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive integer")
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _observation_times(n, spacing, t, rng):
    if t is not None:
        if n is not None or spacing is not None:
            raise ValueError("specify either t or both n and spacing")
        times = _real_array(t, "t")
    else:
        if n is None or spacing is None:
            raise ValueError("specify either t or both n and spacing")
        times = _generate_grid(_positive_int(n, "n"), spacing, rng)
    if times.ndim != 1 or not times.size:
        raise ValueError("observation times must be a nonempty one-dimensional array")
    with np.errstate(over="ignore"):
        gaps = np.diff(times)
    if not (np.isfinite(times).all() and np.isfinite(gaps).all() and (gaps > 0).all()):
        raise ValueError("observation times must be finite and strictly increasing")
    return times


def _generate_grid(n, spacing, rng):
    if not isinstance(spacing, tuple) or len(spacing) != 2:
        raise ValueError("spacing must be a (name, parameter) tuple")
    name, arg = spacing
    if not isinstance(name, str) or name not in ("gamma", "geometric", "single-gap"):
        raise ValueError("spacing name must be 'gamma', 'geometric', or 'single-gap'")
    if n < 2 and name != "geometric":
        raise ValueError(f"{name} spacing requires at least two observations")

    if name == "gamma":
        shape = _finite_real(arg, "Gamma shape")
        if shape <= 0:
            raise ValueError("Gamma shape must be positive")
        try:
            times = varmapack.gamma_times(n, shape=shape, rng=rng)
        except varmapack.VarmapackError as exc:
            raise ValueError(f"could not generate Gamma gaps: {exc}") from exc
        times = 1 + (times / times[-1]) * (n - 1)
        return times

    if name == "geometric":
        miss_prob = _finite_real(arg, "miss_prob")
        if not 0 <= miss_prob < 1:
            raise ValueError("miss_prob must be in [0, 1)")
        if miss_prob == 0 or n == 1:
            return np.arange(1, n + 1, dtype=np.float64)
        # Invert the geometric CDF on {1, 2, ...}.
        with np.errstate(divide="ignore"):
            gaps = 1 + np.floor(np.log1p(-rng.unif(n - 1)) / math.log(miss_prob))
        return np.r_[1.0, 1 + np.cumsum(gaps)]

    if isinstance(arg, tuple) and len(arg) == 2:
        lo, hi = (_positive_int(value, "n_missing bound") for value in arg)
        if lo > hi or hi > np.iinfo(np.int64).max:
            raise ValueError("n_missing bounds must be increasing and fit in int64")
        count = rng.int(a=lo, b=hi) if lo < hi else lo
    else:
        count = _positive_int(arg, "n_missing")
    start = rng.int(a=1, b=n - 1) if n > 2 else 1
    times = np.arange(1, n + 1, dtype=np.float64)
    times[start:] += count
    return times


def _latent_parameters(rho, marginal):
    rho = _finite_real(rho, "rho")
    if not -1 <= rho <= 1:
        raise ValueError("rho must be in [-1, 1]")
    if isinstance(marginal, str) and marginal == "normal":
        return rho, 0.0, 0.0
    if isinstance(marginal, str) and marginal == "lognormal":
        sigma2 = 1.0
        cv = math.sqrt(math.expm1(sigma2))
    elif (
        isinstance(marginal, tuple)
        and len(marginal) == 2
        and isinstance(marginal[0], str)
        and marginal[0] == "lognormal"
        and isinstance(marginal[1], dict)
        and marginal[1].keys() == {"skewness"}
    ):
        g = _finite_real(marginal[1]["skewness"], "skewness")
        if g < 0:
            raise ValueError("skewness must be nonnegative")
        cv = 2 * math.sinh(math.asinh(g / 2) / 3)
        sigma2 = math.log1p(cv * cv)
    else:
        raise ValueError("marginal must be 'normal', 'lognormal', or a lognormal skewness tuple")

    if sigma2 < np.finfo(np.float64).tiny:
        return rho, 0.0, 0.0
    lower = -math.exp(-sigma2)
    if rho < lower:
        raise ValueError(f"rho is infeasible for these lognormal marginals (minimum {lower:g})")
    if rho == lower:
        latent = -1.0
    elif rho == 1:
        latent = 1.0
    elif rho * (cv * cv) < -0.5:
        # Avoid cancellation in 1 + rho * cv**2 near the lower bound.
        latent = math.log(math.exp(-sigma2) + (rho - lower) * (cv * cv)) / sigma2
    else:
        latent = math.log1p(rho * (cv * cv)) / sigma2
    return latent, sigma2, cv


def _make_process(model, tau, drift, theta, lag1, rho):
    cov = [[1.0, rho], [rho, 1.0]]
    if not isinstance(model, str) or model not in ("ou", "ar1", "slepian"):
        raise ValueError("model must be 'ou', 'ar1', or 'slepian'")
    try:
        if model in ("ou", "ar1"):
            if theta is not None or lag1 is not None:
                raise ValueError("theta and lag1 apply only to the Slepian model")
            if drift is not None:
                if tau is not None:
                    raise ValueError("specify either tau or drift, not both")
                drift = _real_array(drift, "drift")
                if drift.shape != (2, 2):
                    raise ValueError("drift must have shape (2, 2)")
                return varmapack.OUModel(drift=drift, cov=cov)
            tau = _real_array((1, 1) if tau is None else tau, "tau")
            if tau.shape != (2,) or not (tau > 0).all():
                raise ValueError("tau must be a pair of positive persistence times")
            return varmapack.OUModel(tau=tau, cov=cov)

        if tau is not None or drift is not None:
            raise ValueError("tau and drift apply only to the OU model")
        if (theta is None) == (lag1 is None):
            raise ValueError("specify exactly one of theta and lag1")
        if lag1 is not None:
            a = _real_array(lag1, "lag1")
            if a.shape not in ((), (2,)) or (np.abs(a) > 0.5).any():
                raise ValueError("lag1 must be a scalar or pair in [-0.5, 0.5]")
            theta = 2 * a / (1 + np.sqrt(1 - 4 * a**2))
        else:
            theta = _real_array(theta, "theta")
        if theta.shape == ():
            theta = np.diag([theta.item(), theta.item()])
        elif theta.shape == (2,):
            theta = np.diag(theta)
        elif theta.shape != (2, 2):
            raise ValueError("theta must be a scalar, pair, or (2, 2) matrix")
        return varmapack.SlepianModel(theta=theta, cov=cov)
    except varmapack.VarmapackError as exc:
        raise ValueError(f"infeasible model parameters or correlation: {exc}") from exc
