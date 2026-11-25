# Project Plan

Work from `~/across`. The project has an unfinished R interface to a Fortran
program that computes bootstrap confidence intervals for correlation between
autocorrelated, cross-correlated time series with potentially uneven time steps.
The intervals are badly biased, especially under the null hypothesis of
independent, uncorrelated series.

Keep local Codex plans in the Git-ignored `codex-plans/` directory.
The separate `studies/` repository owns the calibration analyses,
`studies/simul/` runner, and `studies/archive/2026-summer/` archive.
The runner is synced to `~/across/simul` on Elja; remote paths are unchanged.

## Objectives

a. Implement a bias-correction backend.
b. Turn the program into an R package for CRAN.
c. Create a Python package for PyPI.
d. Write an article about the method.

## Bias-Correction Dimensions

- Sample size (`n`).
- Persistence (`tau_x` and `tau_y`).
- Time-step variability (Gamma spacing).
- Marginal distribution (normal and log-normal).
- True correlation (`rho`).
- VARMA kind (`p=1, q=0`, and `p=q=2`).
