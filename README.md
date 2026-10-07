# aersn: affine-equivariant adjusted-range self-normalization

**Version 0.2.3**

`aersn` provides tests, confidence intervals and joint confidence regions
for parameters estimated from dependent time series. Start with a sample
mean, a fitted regression, or your own estimate and influence contributions.
Use the same fitted object to compare adjusted-range inference with LDL
partial prewhitening, Shao self-normalization, heteroskedasticity and
autocorrelation consistent (HAC) covariance estimation, Bartlett fixed-b,
and equal-weighted cosine (EWC) inference.

The main method, affine-equivariant adjusted-range self-normalization,
constructs one joint region from the ranges of the centered influence path
in all directions. Changing the units or applying a nonsingular linear
transformation changes the region in the corresponding way. It needs no
kernel or bandwidth and reduces to scalar adjusted-range inference when
there is one parameter. The construction uses the convex hull of path
increments; tests and simultaneous intervals use its gauge and support
function. See [Hong, Lin, Linton, Newey and Sun (2026), Cambridge Working Papers
in Economics No. 2678](https://www.janeway.econ.cam.ac.uk/publication/affine-equivariant-adjusted-range-self-normalization)
and Hong, Linton, McCabe, Sun and Wang (2024, *Journal of Econometrics*).

A supplied or consistently estimated common variance-accumulation profile
can be used with the hull, LDL and Shao methods. For Shao, the profile option also
changes the integration weights; the adjusted-range construction uses path
ranges. These options require the common-profile conditions described in
`help("aersn_profile")`.

## Installation

For a CRAN release, install with `install.packages("aersn")`.
For a source archive supplied by the authors, use the following commands.
Only `lpSolve` and `sandwich` are needed beyond the packages that ship with R.

```r
install.packages(c("lpSolve", "sandwich"))
install.packages("/path/to/aersn_0.2.3.tar.gz", repos = NULL, type = "source")
```

The source tarball includes built tutorials. After installation, open them
with `vignette(package = "aersn")`; no vignette-building step is needed.
Only developers rebuilding tutorials from the repository need `knitr` and
`rmarkdown`.

## A minimal example

```r
library(aersn)
set.seed(1)
n <- 300
e <- matrix(rnorm(2 * n), n, 2)
Y <- e
for (t in 2:n) Y[t, ] <- 0.5 * Y[t - 1, ] + e[t, ]   # bivariate AR(1)

fit <- aersn_mean(Y, names = c("m1", "m2"))          # psi_t = Y_t - Ybar
aersn_test(fit, null = c(0, 0))                      # increment-hull test
confint(fit)                                         # simultaneous intervals
aersn_contrast(fit, c(1, -1))                        # a linear contrast
plot(aersn_region(fit))                              # joint region
```

## Six methods on one estimate

| `method` | Normalizer | Tuning | Reference law |
|---|---|---|---|
| `"hull"` (default) | increment hull of the centered path | none | simulated Brownian gauge law |
| `"ldl"` | componentwise adjusted ranges after lag-zero prewhitening | coordinate order | simulated independent-component law |
| `"shao"` | integrated outer product of the path | integration rule | simulated Brownian quadratic law |
| `"hac"` | kernel long-run covariance estimate | kernel, bandwidth | chi-squared |
| `"fixedb"` | Bartlett estimate with bandwidth fraction *b* | *b* | simulated fixed-b law |
| `"ewc"` | equal-weighted cosine estimate | number of terms | scaled F |

```r
aersn_compare(fit, null = c(0, 0))

aersn_test(fit, method = "shao")
aersn_test(fit, method = "hac", kernel = "Parzen", bandwidth = "andrews")
aersn_test(fit, method = "fixedb", b = 0.5)
aersn_test(fit, method = "ewc", nu = 20)
confint(fit, method = "ldl")
```

Each result records the tuning values actually used, including bandwidths
selected from the data and fixed-b fractions rounded to an integer lag
bandwidth. Statistics and critical values are on different scales across
methods and are not comparable as numbers; p-values, decisions and interval
widths are. The vignette *Comparing inference methods on one estimate* works
through this.

## Estimators other than the mean

Supply an estimate and its observation-level influence contributions to
`aersn()`, or use a model interface that computes them: `aersn_lm()` for
least squares, `aersn_gmm()` for smooth generalized method of moments
including instrumental variables, and `aersn_mle()` for conditional
likelihood scores. For a parameter that is a function of a larger estimated
vector, `aersn_target()` applies the Jacobian, so that jointly estimated
nuisance coefficients keep their first-order effect. All interfaces feed the
same core construction.

## Documentation

Seven vignettes: scalar mean inference; multivariate mean inference and
linear contrasts; supplied influence contributions; linear regression and
smooth generalized method of moments; conditional likelihood scores and
variance-accumulation profiles; reference distributions and reproducibility;
and comparing inference methods. Use `help(package = "aersn")` for the
function index and `vignette(package = "aersn")` for the installed tutorials.

## What the methods assume

All six need the influence contributions to satisfy a functional central
limit theorem with a nonsingular long-run covariance matrix, and the
estimator to be asymptotically linear in them. The package cannot verify
those conditions; each method's help page states what else it needs. The
componentwise method's reference law requires the transformed long-run
covariance to be diagonal in the limit, which diagonalizing the sample
lag-zero covariance does not deliver, and its statistic is not affine
equivariant. HAC inference relies on the conditions under which the estimate
is consistent. Fixed-b and EWC reference laws are fixed-smoothing asymptotic
laws, not exact finite-sample distributions. Strong persistence can cause
finite-sample size distortion; a matched reference grid does not remove it.

## Scope

Not included: the Kolmogorov-Smirnov type structural-break test of Hong et
al. (2024) and the autocorrelation tests of Sun, Zhu and Linton (2025). The
comparator methods here are not adjusted-range versions of those procedures.

## License

MIT, held by the five authors. Maintainer: Jiajing Sun
(jiajing.sun@gmail.com).
