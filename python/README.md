# autocross for Python

Install from the repository root with a C compiler and GNU Fortran on
your PATH:

```sh
python -m pip install ./python
```

Generate a stationary pair on a Gamma-spaced observation grid:

```python
from autocross import testcase

t, x, y = testcase(100, model="ou", spacing=("gamma", 16), tau=(1, 2))
```

Compute the uncorrected PearsonT3 confidence interval:

```python
from autocross import pearson_t3

result = pearson_t3(t, x, y)
print(result.r, result.ci)
```

See `help(testcase)` and `help(pearson_t3)` for the full interfaces. The
Python package calls Fortran directly and does not require R or a separate
native Randompack installation. For reproducible calls, pass a seeded
`randompack.Rng`. Use separate processes to parallelize simulations.

To run the tests and doctests: `python -m pip install './python[test]'`, then
`python -m pytest --import-mode=importlib python/tests --pyargs autocross --doctest-modules`.
The optional standalone comparison uses `AUTOCROSS_TEST_FORTRAN_DRIVER`
pointing to an executable built from `tests/pearson_t3_reference.f90`, the
legacy Fortran sources, and the native Randompack Fortran module and library.
