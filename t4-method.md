# PearsonT4 Bootstrap Calibration

Let the observed paired time series be

$$
D=\{(t_i,x_i,y_i): i=1,\ldots,n\},
$$

and let the ordinary Pearson correlation be

$$
r=\operatorname{cor}(x,y), \qquad z=\operatorname{atanh}(r).
$$

The point estimate is the usual Pearson sample correlation. The bootstrap machinery
is used only to construct the confidence interval.

## Notation

We use uppercase symbols for bootstrap counts. For T3, the usual nested
bootstrap counts are

$$
B_1,\quad B_2.
$$

For T4, the null-calibration counts are

$$
B_1,\quad B_2,\quad B_3,
$$

where \(B_1\) is the number of simulated null datasets, \(B_2\) is the number
of first-level bootstrap samples within each null dataset, and \(B_3\) is the
number of second-level bootstrap samples within each first-level sample.

## T3 Bootstrap

PearsonT3 uses a nested moving-block bootstrap. Let \(B_1\) be the number of
first-level bootstrap samples and \(B_2\) the number of second-level bootstrap
samples.

For \(j=1,\ldots,B_1\), draw a first-level bootstrap sample

$$
D_j^*
$$

from the observed data and compute

$$
r_j^*=\operatorname{cor}(x_j^*,y_j^*), \qquad
z_j^*=\operatorname{atanh}(r_j^*).
$$

For each first-level sample \(D_j^*\), draw second-level bootstrap samples

$$
D_{jk}^{**}, \qquad k=1,\ldots,B_2,
$$

and compute

$$
r_{jk}^{**}=\operatorname{cor}(x_{jk}^{**},y_{jk}^{**}), \qquad
z_{jk}^{**}=\operatorname{atanh}(r_{jk}^{**}).
$$

The conditional standard error associated with first-level sample \(j\) is estimated
as

$$
s_j = \operatorname{sd}_{k=1}^{B_2}(z_{jk}^{**}).
$$

The T3 interval is a Studentized bootstrap interval on Fisher scale, transformed
back to correlation scale:

$$
\operatorname{CI}_{T3}
= \tanh\left(z \pm \lambda_{T3}s_{\mathrm{obs}}\right),
$$

where \(s_{\mathrm{obs}}\) is the estimated standard error for the observed data and
\(\lambda_{T3}\) is obtained from the T3 bootstrap calibration.

## T4 Null Calibration

PearsonT4 keeps the Pearson point estimate, but replaces the T3 confidence
interval construction by a null-calibrated interval. The goal is that a nominal
95 % interval has approximately 5 % two-sided rejection probability when the
true correlation is zero.

Let \(B_1\) be the number of null calibration datasets. For
\(m=1,\ldots,B_1\), generate an independent-null dataset

$$
D_m^0
$$

by drawing moving-block bootstrap samples from \(x\) and \(y\) independently.
This preserves the marginal data values and serial dependence approximately,
while breaking the contemporaneous association so that the null correlation is

$$
\rho=0.
$$

Concretely, one null sample \(D^0\) is constructed from the observed paired
series

$$
D=\{(x_i,y_i)\}_{i=1}^n
$$

as follows. First draw a moving-block bootstrap sample from the observed \(x\)
series only,

$$
x^0=(x_1^0,\ldots,x_n^0).
$$

Then independently draw a moving-block bootstrap sample from the observed \(y\)
series only,

$$
y^0=(y_1^0,\ldots,y_n^0).
$$

Finally pair those two independently resampled series by index:

$$
D^0=\{(x_i^0,y_i^0)\}_{i=1}^n.
$$

Thus a null sample \(D^0\) is not generated from an AR model. It is made
directly from the observed data by independently block-resampling \(x\) and
\(y\). The independent block draws impose the null hypothesis \(\rho=0\), while
the block structure attempts to preserve each series' own autocorrelation.

For each null dataset \(D_m^0\), run a nested bootstrap with calibration counts
\(B_2\) and \(B_3\):

$$
D_m^0
\longrightarrow
D_{mj}^{0*}, \qquad j=1,\ldots,B_2,
$$

and

$$
D_{mj}^{0*}
\longrightarrow
D_{mjk}^{0**}, \qquad k=1,\ldots,B_3.
$$

For each second-level null bootstrap sample, compute

$$
r_{mjk}^{0**}
= \operatorname{cor}(x_{mjk}^{0**},y_{mjk}^{0**}),
\qquad
z_{mjk}^{0**}
= \operatorname{atanh}(r_{mjk}^{0**}).
$$

For each first-level null sample, estimate

$$
s_{mj}^0
=
\operatorname{sd}_{k=1}^{B_3}(z_{mjk}^{0**}).
$$

The T4 implementation summarizes these first-level standard errors into a
dataset-level null standard error,

$$
s_m^0
=
\frac{1}{B_2}
\sum_{j=1}^{B_2} s_{mj}^0.
$$

It then summarizes the null-dataset standard errors into a single T4 scale,

$$
s_{T4}
=
\frac{1}{B_1}
\sum_{m=1}^{B_1} s_m^0.
$$

Let

$$
r_m^0=\operatorname{cor}(x_m^0,y_m^0),
\qquad
z_m^0=\operatorname{atanh}(r_m^0).
$$

For candidate lower and upper endpoint values \(\lambda_L\) and \(\lambda_U\),
the null interval for dataset \(m\) is

$$
\operatorname{CI}_m^0(\lambda_L,\lambda_U)
=
\left[
\tanh\left(z_m^0 + t(\lambda_L)s_{T4}\right),
\tanh\left(z_m^0 - t(\lambda_U)s_{T4}\right)
\right],
$$

where \(t(\lambda)\) is the Student-\(t\) quantile used by PearsonT. The two
null tail error rates are estimated separately:

$$
\widehat{\alpha}_{L}(\lambda_L)
=
\frac{1}{B_1}
\sum_{m=1}^{B_1}
I\{L_m^0(\lambda_L)>0\},
$$

and

$$
\widehat{\alpha}_{U}(\lambda_U)
=
\frac{1}{B_1}
\sum_{m=1}^{B_1}
I\{U_m^0(\lambda_U)<0\}.
$$

For a 95 % two-sided confidence interval, PearsonT4 uses the one-tail level

$$
\alpha = 0.025.
$$

It chooses \(\lambda_L\) and \(\lambda_U\) so that

$$
\widehat{\alpha}_{L}(\lambda_L) \approx \alpha,
\qquad
\widehat{\alpha}_{U}(\lambda_U) \approx \alpha.
$$

The final T4 interval for the observed data is

$$
\operatorname{CI}_{T4}
=
\left[
\tanh\left(z + t(\lambda_L)s_{T4}\right),
\tanh\left(z - t(\lambda_U)s_{T4}\right)
\right].
$$

Thus the reported T4 interval uses the same null-path scale \(s_{T4}\) that was
used to calibrate the null rejection probability. It does not use the observed-data
T3 bootstrap standard error. Because \(\lambda_L\) and \(\lambda_U\) are chosen
separately, the T4 interval need not be symmetric on Fisher scale.

## Computational Cost

The null-calibration part has cost approximately

$$
B_1 B_2 B_3.
$$

For example, with

$$
B_1=1000,\qquad B_2=100,\qquad B_3=50,
$$

the null calibration costs roughly

$$
1000\cdot 100\cdot 50 = 5\times 10^6
$$

correlation calculations.
