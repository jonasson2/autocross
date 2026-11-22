from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import os
from statistics import NormalDist
import subprocess
import warnings

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest
from randompack import Rng

from autocross import (
    AutocrossError, CalibrationWarning, DegenerateBootstrapWarning, PersistenceWarning,
    T3Result, T3Settings, pearson_t3, testcase as make_case,
)
from autocross._t3 import _pearson_t3_levels


def seeded(seed=123):
    rng = Rng("x256++simd")
    rng.seed(seed)
    return rng


@pytest.fixture
def data():
    return make_case(80, model="ou", spacing=("gamma", 16), rho=0.3, rng=seeded())


def test_multiple_confidences_reuse_bootstrap(data):
    settings = T3Settings(b1=100, b2=60, n_lambda=101, block_length=4)
    levels = (0.95, 0.64, 0.98, 0.9)
    stream = seeded()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CalibrationWarning)
        results = _pearson_t3_levels(*data, levels=levels, rng=stream, settings=settings, diagnostics=True)
        for level, result in zip(levels, results):
            reference_rng = seeded()
            reference = pearson_t3(*data, confidence=level, rng=reference_rng, settings=settings, diagnostics=True)
            assert result.confidence == level
            assert_allclose(result.ci, reference.ci, rtol=0, atol=2e-15)
            assert result.diagnostics.selected_lambda == reference.diagnostics.selected_lambda
            assert result.diagnostics.calibration_curve is results[0].diagnostics.calibration_curve
    assert_array_equal(stream.unif(20), reference_rng.unif(20))


def test_repeated_calls_and_diagnostics(data, capfd, monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    settings = T3Settings(b1=100, b2=60, n_lambda=101, block_length=4)
    original = pearson_t3(*data, rng=seeded(), settings=settings, diagnostics=True)
    assert isinstance(original, T3Result)
    assert original.block_length == 4
    assert original.ci[0] <= original.r <= original.ci[1]
    assert original.taux > 0 and original.tauy > 0
    assert original.settings is settings
    assert original.confidence == 0.95
    diag = original.diagnostics
    assert_allclose(diag.se_z, diag.bootstrap_z.std(ddof=1))
    assert diag.selected_lambda == diag.lambda_grid[np.argmin(abs(diag.calibration_curve - 0.95))]
    assert diag.lambda_grid.shape == diag.calibration_curve.shape == (101,)
    assert diag.bootstrap_z.shape == diag.inner_se_z.shape == (100,)
    assert not diag.bootstrap_z.flags.writeable
    snapshot = diag.bootstrap_z.copy()
    for key in ("P3_SEED", "P3_NO_CALIBRATION", "P3_NO_FISHER", "P3_ORIG_PEARSN",
                "P3_L_MBB", "P3_POLROM", "P4_POLROM", "P3_PROFILE", "P3_BOOT_SLOW"):
        monkeypatch.setenv(key, "1")
    pearson_t3(*(a[:70] for a in data), rng=seeded(99),
               settings=T3Settings(b1=80, b2=40, detrend="linear"))
    repeat = pearson_t3(*data, rng=seeded(), settings=settings)
    assert repeat.diagnostics is None
    assert (repeat.r, repeat.ci, repeat.taux, repeat.tauy) == (
        original.r, original.ci, original.taux, original.tauy)
    assert_array_equal(diag.bootstrap_z, snapshot)
    assert capfd.readouterr() == ("", "")
    assert not list(tmp_path.iterdir())


def test_strided_inputs_and_rng(data):
    t, x, y = (a[::2] for a in data)
    saved = [a.copy() for a in (t, x, y)]
    settings = T3Settings(b1=100, b2=50, block_length=3)
    rng = seeded()
    with pytest.warns(PersistenceWarning):
        result = pearson_t3(t, x, y, rng=rng, settings=settings)
        again = pearson_t3(t, x, y, rng=rng, settings=settings)
    assert result.ci != again.ci
    ref = seeded()
    ref.int(size=settings.b1 * (settings.b2 + 1) * 14, a=1, b=38, dtype=np.int32)
    ref.int(size=settings.b1 * (settings.b2 + 1) * 14, a=1, b=38, dtype=np.int32)
    assert_array_equal(rng.int(size=30, a=1, b=38), ref.int(size=30, a=1, b=38))
    for a, b in zip((t, x, y), saved):
        assert_array_equal(a, b)


def test_linear_detrending(data):
    t, x, y = data
    settings = T3Settings(b1=100, b2=60, detrend="linear", block_length=3)
    result = pearson_t3(t, x + 3*t, y - 2*t, settings=settings, rng=seeded())
    baseline = pearson_t3(t, x, y, settings=settings, rng=seeded())
    assert_allclose(result.r, baseline.r, atol=1e-13)
    assert_allclose(result.ci, baseline.ci, atol=1e-13)
    with pytest.raises(ValueError, match="constant after detrending"):
        pearson_t3(np.arange(x.size), np.arange(x.size), y, settings=settings)


def test_student_quantile_for_long_series():
    n = 500
    t, x, y = make_case(model="ou", t=range(n), tau=(3, 6), rho=0.3, rng=seeded())
    settings = T3Settings(b1=100, b2=60, n_lambda=40, block_length=10)
    result = pearson_t3(t, x, y, confidence=0.85, settings=settings, rng=seeded(), diagnostics=True)
    d = 2*n - 5
    z = NormalDist().inv_cdf(result.diagnostics.selected_lambda)
    quantile = z + (z**3 + z)/(4*d) + (5*z**5 + 16*z**3 + 3*z)/(96*d**2)
    quantile += (3*z**7 + 19*z**5 + 17*z**3 - 15*z)/(384*d**3)
    quantile += (79*z**9 + 776*z**7 + 1482*z**5 - 1920*z**3 - 945*z)/(92160*d**4)
    expected = np.tanh(np.arctanh(result.r) + np.array([1, -1])*quantile*result.diagnostics.se_z)
    assert_allclose(result.ci, expected, atol=1e-8, rtol=0)


def test_calibration_warning(data):
    settings = T3Settings(b1=100, b2=30, n_lambda=1, block_length=3)
    with pytest.warns(CalibrationWarning, match="outside the calibration coverage range"):
        result = pearson_t3(*data, settings=settings, rng=seeded())
    assert_allclose(result.ci, [result.r, result.r])
    with warnings.catch_warnings():
        warnings.simplefilter("error", CalibrationWarning)
        with pytest.raises(CalibrationWarning):
            pearson_t3(*data, settings=settings, rng=seeded(), diagnostics=True)


def test_disputed_persistence_refit_matches_bounded_least_squares():
    # Rounded uneven-grid example from the n=20 null pilot. Only x needs a refit.
    t = np.array([1, 2.599, 3.184, 4.176, 5.664, 6.665, 7.805, 8.417, 9.468, 10.665,
                  11.691, 12.917, 14.004, 15.047, 16.029, 16.815, 17.683, 18.461, 19.475, 20])
    x = np.array([1.312, -.593, -.601, -.950, .541, .360, -.280, .222, -.368, -.297,
                  .315, -.760, -.858, .225, .407, 1.106, -.006, -1.160, -.424, .890])
    v = (x - x.mean())[::-1]
    d = np.diff(t)[::-1]
    q = np.linspace(0, 1, 20001)[:, None]
    objective = ((v[1:] - v[:-1] * q**(d / d.mean()))**2).sum(axis=1)
    q_raw = q[objective.argmin(), 0]
    q_corrected = ((len(x) - 1) * q_raw + 1) / (len(x) - 4)
    for scale in (1, 10):
        with pytest.warns(PersistenceWarning, match="bounded persistence refit for x$"):
            result = pearson_t3(scale*t, x, np.sin(t), settings=T3Settings(b1=128, b2=64),
                                diagnostics=True, rng=seeded())
        assert result.diagnostics.persistence_refit == (True, False)
        assert abs(np.exp(-scale*d.mean()/result.taux) - q_corrected) < 5e-5
        assert np.isfinite(result.ci).all()


def test_persistence_and_degenerate_warnings():
    t = np.arange(20)
    x = (-1.0)**t
    y = np.sin(t * 0.43)
    with pytest.warns(PersistenceWarning, match="for x"):
        pearson_t3(t, x, y, rng=seeded(), settings=T3Settings(b1=100, b2=30))
    x = np.r_[np.zeros(19), 1.]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pearson_t3(t, x, y, rng=seeded(), settings=T3Settings(b1=100, b2=30, block_length=1))
    assert sum(type(w.message) is DegenerateBootstrapWarning for w in caught) == 1


def test_failure_does_not_poison_next_call(data):
    t, x, y = data
    with pytest.raises(AutocrossError, match="could not estimate persistence"):
        pearson_t3(np.arange(t.size)*1e-310, x, y, rng=seeded())
    assert np.isfinite(pearson_t3(*data, rng=seeded(), settings=T3Settings(b1=100, b2=60)).r)


def test_process_parallel_calls(data):
    settings = T3Settings(b1=128, b2=64)
    with ProcessPoolExecutor(2, mp_context=get_context("spawn")) as pool:
        calls = [pool.submit(pearson_t3, *data, settings=settings) for _ in range(2)]
        results = [call.result(timeout=30) for call in calls]
    assert results[0].r == results[1].r
    assert all(result.ci[0] <= result.r <= result.ci[1] for result in results)


@pytest.mark.parametrize("kwargs, error", [
    ({"settings": {}}, TypeError), ({"confidence": "0.9"}, TypeError),
    ({"confidence": 1}, ValueError), ({"rng": 123}, TypeError),
    ({"diagnostics": 1}, TypeError), ({"settings": T3Settings(block_length=80)}, ValueError),
])
def test_invalid_options(data, kwargs, error):
    with pytest.raises(error):
        pearson_t3(*data, **kwargs)


def test_invalid_data(data):
    t, x, y = data
    for args in ((t, x[:-1], y), (t[:9], x[:9], y[:9]), (t[::-1], x, y),
                 (t, np.zeros_like(x), y), (t, x*float("nan"), y)):
        with pytest.raises(ValueError):
            pearson_t3(*args)


@pytest.mark.parametrize("key, value, error", [
    ("b1", True, TypeError), ("b2", 1, ValueError), ("n_lambda", 0, ValueError),
    ("block_length", 0, ValueError), ("detrend", "polynomial", ValueError),
])
def test_invalid_settings(key, value, error):
    with pytest.raises(error):
        T3Settings(**{key: value})


@pytest.mark.parametrize("even, linear, block", [(False, False, None), (True, True, 3)])
def test_legacy_fortran_driver(even, linear, block, data):
    driver = os.environ.get("AUTOCROSS_TEST_FORTRAN_DRIVER")
    if not driver:
        pytest.skip("standalone Fortran reference driver not configured")
    t, x, y = data
    if even:
        t = np.arange(t.size, dtype=float)
    settings = T3Settings(b1=128, b2=64, n_lambda=101, block_length=block,
                          detrend="linear" if linear else "mean")
    env = {k: v for k, v in os.environ.items() if not k.startswith(("P3_", "P4_"))}
    env["P3_SEED"] = "123"
    if block is not None:
        env["P3_L_MBB"] = str(block)
    header = f"{len(t)} 128 64 101 0.95 {'l' if linear else 'm'}\n"
    rows = "\n".join(" ".join(format(v, ".17g") for v in row) for row in zip(t, x, y))
    stdout = subprocess.run([driver], input=header+rows+"\n", text=True, env=env,
                            capture_output=True, check=True).stdout
    expected = next(line for line in stdout.splitlines() if line.startswith("RESULT"))
    result = pearson_t3(t, x, y, settings=settings, rng=seeded())
    assert_allclose([result.r, *result.ci, result.taux, result.tauy, result.block_length],
                    [float(v) for v in expected.split()[1:]], atol=1e-12, rtol=1e-12)
