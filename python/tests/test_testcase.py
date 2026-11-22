import math

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from randompack import Rng
from varmapack import OUModel, SlepianModel, gamma_times

import autocross


def seeded():
    rng = Rng()
    rng.seed(123)
    return rng


def test_explicit_grid_and_defaults():
    t, x, y = autocross.testcase(model="ou", t=range(2, 22, 2), rng=seeded())
    assert_array_equal(t, range(2, 22, 2))
    assert all(a.shape == (10,) and a.dtype == np.float64 for a in (t, x, y))
    expected = OUModel(tau=[1, 1], cov=np.eye(2)).sim(t, rng=seeded())[0]
    assert_array_equal(np.column_stack((x, y)), expected)


def test_gamma_grid_and_alias():
    settings = dict(spacing=("gamma", 16), tau=(1, 2), rho=0.3)
    actual = autocross.testcase(100, model="ar1", rng=seeded(), **settings)
    for a, b in zip(actual, autocross.testcase(100, model="ou", rng=seeded(), **settings)):
        assert_array_equal(a, b)
    t, x, y = actual
    assert t[0] == 1 and t[-1] == 100
    assert (np.diff(t) > 0).all()
    rng = seeded()
    raw = gamma_times(100, shape=16, rng=rng)
    assert_allclose(t, 1 + raw * (99 / raw[-1]))
    expected = OUModel(tau=[1, 2], cov=[[1, 0.3], [0.3, 1]]).sim(t, rng=rng)[0]
    assert_array_equal(np.column_stack((x, y)), expected)


def test_rng_advances():
    rng = seeded()
    first = autocross.testcase(model="ou", t=[1, 2], rng=rng)
    second = autocross.testcase(model="ou", t=[1, 2], rng=rng)
    assert not np.array_equal(first[1], second[1])


def test_geometric_grid():
    t, *_ = autocross.testcase(5000, model="ou", spacing=("geometric", 0.25), rng=seeded())
    gaps = np.diff(t)
    assert t[0] == 1
    assert (gaps >= 1).all() and (gaps == np.floor(gaps)).all()
    assert gaps.mean() == pytest.approx(1 / 0.75, abs=0.025)
    t, *_ = autocross.testcase(5, model="ou", spacing=("geometric", 0), rng=seeded())
    assert_array_equal(t, np.arange(1, 6))


@pytest.mark.parametrize("missing", [5, (3, 8), (4, 4)])
def test_single_gap(missing):
    t, *_ = autocross.testcase(50, model="ou", spacing=("single-gap", missing), rng=seeded())
    gaps = np.diff(t)
    assert t[0] == 1 and t.size == 50
    assert (gaps == 1).sum() == 48
    skipped = gaps.max() - 1
    if isinstance(missing, tuple):
        assert missing[0] <= skipped <= missing[1]
    else:
        assert skipped == missing


def test_coupled_ou():
    t = np.array([1.0, 1.25, 5.0])
    drift = [[1.0, 0.3], [0.1, 0.5]]
    _, x, y = autocross.testcase(model="ar1", t=t, drift=drift, rng=seeded())
    expected = OUModel(drift=drift, cov=np.eye(2)).sim(t, rng=seeded())[0]
    assert_array_equal(np.column_stack((x, y)), expected)


@pytest.mark.parametrize("params, coefficient", [
    ({"theta": 0.5}, np.diag([0.5, 0.5])),
    ({"theta": (0.5, 0.25)}, np.diag([0.5, 0.25])),
    ({"lag1": (0.4, -0.5)}, np.diag([0.5, -1.0])),
    ({"theta": [[0.5, 0.1], [0.0, 0.25]]}, [[0.5, 0.1], [0.0, 0.25]]),
])
def test_slepian(params, coefficient):
    t = np.array([1.0, 1.5, 4.5, 5.0])
    _, x, y = autocross.testcase(model="slepian", t=t, rng=seeded(), **params)
    process = SlepianModel(theta=coefficient, cov=np.eye(2))
    assert_allclose(np.column_stack((x, y)), process.sim(t, rng=seeded())[0])


def test_lognormal_default_sd_and_gaussian_limit():
    kwargs = dict(model="ou", t=range(20), rho=0.3)
    normal = autocross.testcase(rng=seeded(), **kwargs)
    zero = autocross.testcase(rng=seeded(), marginal=("lognormal", {"skewness": 0}), **kwargs)
    small = autocross.testcase(rng=seeded(), marginal=("lognormal", {"skewness": 1e-4}), **kwargs)
    for a, b, c in zip(normal, zero, small):
        assert_array_equal(a, b)
        assert_allclose(a, c, atol=2e-4)

    _, x, y = autocross.testcase(rng=seeded(), marginal="lognormal", **kwargs)
    cv2 = math.expm1(1)
    r = math.log1p(kwargs["rho"] * cv2)
    cov = [[1, r], [r, 1]]
    z = OUModel(tau=[1, 1], cov=cov).sim(normal[0], rng=seeded())[0]
    assert_allclose(np.column_stack((x, y)), (np.exp(z - 0.5) - 1) / math.sqrt(cv2))


@pytest.mark.parametrize("rho", [0.4, -0.4])
def test_lognormal_population_targets(rho):
    _, x, y = autocross.testcase(
        model="ou", t=range(100_000), tau=(0.1, 0.1), rho=rho, rng=seeded(),
        marginal=("lognormal", {"skewness": 2}),
    )
    assert np.corrcoef(x, y)[0, 1] == pytest.approx(rho, abs=0.015)
    for a in (x, y):
        assert a.mean() == pytest.approx(0, abs=0.02)
        assert a.std() == pytest.approx(1, abs=0.03)
        assert np.mean(((a - a.mean()) / a.std())**3) == pytest.approx(2, abs=0.15)


@pytest.mark.parametrize("kwargs, message", [
    ({"t": [1, 2], "n": 2}, "either t"),
    ({"n": 2}, "either t"),
    ({"spacing": ("gamma", 16)}, "either t"),
    ({"n": True, "spacing": ("gamma", 16)}, "positive integer"),
    ({"n": 1, "spacing": ("gamma", 16)}, "at least two"),
    ({"n": 2, "spacing": "gamma"}, "tuple"),
    ({"n": 2, "spacing": ("gamma", 0)}, "positive"),
    ({"n": 2, "spacing": ("geometric", 1)}, "miss_prob"),
    ({"n": 2, "spacing": ("single-gap", (8, 3))}, "bounds"),
    ({"t": [[1, 2]]}, "one-dimensional"),
    ({"t": [1, 1]}, "strictly increasing"),
    ({"t": [1, np.inf]}, "finite real"),
    ({"t": [1, 2j]}, "finite real"),
    ({"tau": (1, 0)}, "positive persistence"),
    ({"tau": (1, 1), "drift": [[1, 0], [0, 1]]}, "either tau or drift"),
    ({"drift": [[-1, 0], [0, 1]]}, "infeasible model"),
    ({"theta": 0.5}, "only to the Slepian"),
    ({"marginal": ("lognormal", 2)}, "skewness tuple"),
    ({"marginal": ("lognormal", {"skewness": -1})}, "nonnegative"),
    ({"marginal": "lognormal", "rho": -0.5}, "rho is infeasible"),
    ({"tau": (1, 10), "rho": 0.9}, "infeasible model"),
])
def test_invalid_ou_inputs(kwargs, message):
    if "n" not in kwargs and "spacing" not in kwargs:
        kwargs = {"t": [1, 2], **kwargs}
    with pytest.raises(ValueError, match=message):
        autocross.testcase(model="ou", rng=seeded(), **kwargs)


def test_invalid_slepian_inputs():
    with pytest.raises(ValueError, match="exactly one"):
        autocross.testcase(model="slepian", t=[1, 2])
    with pytest.raises(ValueError, match="scalar or pair"):
        autocross.testcase(model="slepian", t=[1, 2], lag1=0.6)
    with pytest.raises(ValueError, match="exactly one"):
        autocross.testcase(model="slepian", t=[1, 2], theta=0.5, lag1=0.4)
