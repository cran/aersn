# aersn 0.2.3

* Use pivoted QR preconditioning for the gauge linear program to retain
  accuracy for nearly collinear, full-rank influence paths. Check the
  returned direction for feasibility and clarify solver failure messages.
* Standardize the GMM first-order-condition warning by moment variation;
  record the residual and numerical tolerance without treating the check
  as a test of asymptotic assumptions.
* Clarify the registry's MIT license, installation instructions, scalar
  reference memory costs, and the target dimension used by HAC adjustment.
* Repair working-paper links, author initials and PDF guide output wrapping.
* Add regression tests for ill-conditioning, pivot/direction mapping,
  test-region consistency and optimizer-based GMM estimates.

# aersn 0.2.2

* Update the main-method citation to Cambridge Working Papers in Economics
  No. 2678 and Janeway Institute Working Paper No. 2637, with the public URL.
* Prepare the first CRAN submission; no changes to numerical methods or APIs.

# aersn 0.2.1

* Bind reusable normalizers to their influence contributions and profile;
  reject mismatched data and tuning supplied with a prebuilt object.
* Apply calendar-time restrictions to direct HAC, fixed-b and EWC
  constructors as well as the common interface.
* Check quadratic joint regions and intervals for singularity, assessing
  dependence after coordinate scaling so changes of units remain usable.
* Correct the returned LDL transformation for a requested coordinate order,
  reject fractional permutations, and improve LDL scaling checks.
* Validate fixed-b reference bandwidth/grid consistency and reference arguments.
* Report separately selected tuning for each marginal contrast.
* Stabilize quadratic spectral weights near zero; clarify bandwidth versus lag,
  the asymptotic EWC reference, plot labels and installation instructions.

# aersn 0.2.0

Local review candidate. Not released, not submitted.

## Five further inference methods

Inference is no longer restricted to the affine-equivariant adjusted-range
increment hull, which remains the default and is unchanged. Five methods can
now be applied to the same fit, the same influence contributions, the same
null value and the same level, through a `method` argument on
`aersn_test()`, `confint()`, `aersn_contrast()`, `aersn_region()` and
`plot()`, and side by side through the new `aersn_compare()`:

* `"ldl"`, componentwise adjusted ranges after lag-zero partial prewhitening
  by an LDL factorization (`aersn_ldl_normalizer()`, `aersn_ldl()`);
* `"shao"`, quadratic self-normalization, with calendar-time or
  variance-profile integration (`aersn_shao_normalizer()`);
* `"hac"`, kernel long-run covariance estimation with the Bartlett, Parzen
  and quadratic spectral kernels, four bandwidth rules, explicit bandwidths
  or lag truncations, optional VAR(1) prewhitening and an optional
  finite-sample multiplier (`aersn_hac_lrv()`);
* `"fixedb"`, Bartlett fixed-b (`aersn_fixed_b_normalizer()`);
* `"ewc"`, the equal-weighted cosine method (`aersn_ewc_lrv()`).

`aersn_normalizer()` builds the normalizer of any method and reports the
tuning values actually used, including bandwidths selected from the data and
fixed-b fractions rounded to an integer lag truncation.

## Reference laws

`aersn_reference()` gains `statistic` and `args` arguments and simulates the
matched-grid law of the componentwise, quadratic and fixed-b statistics as
well as the increment-hull gauge. `aersn_parametric_reference()` provides the
chi-squared law used with HAC and the scaled F law used with EWC. Cache keys,
and the check that a statistic is combined with its own law, now include the
method and every tuning value that changes the distribution, so a law
simulated for one method or one bandwidth fraction cannot be used with
another. The scalar componentwise law reuses the scalar hull draws, so the
two agree exactly.

## Verified identities

The test suite checks, among others, that the Bartlett fixed-b estimate at
`b = 1` is exactly twice the quadratic self-normalizer and that the two tests
then agree on common reference draws; that the HAC estimate and its Wald
statistic reproduce `sandwich::vcovHAC()` to machine precision for all three
kernels, with and without prewhitening; that the Andrews bandwidths for the
Parzen and quadratic spectral kernels stand in the ratio 2.6614/1.3221; that
fixed-b at a realized fraction `m/n` equals the Bartlett HAC estimate at lag
truncation `m - 1`, the two conventions differing by one; and that the
componentwise statistic is unchanged by the choice of divisor in the lag-zero
covariance and by the use of the Cholesky factor in place of the unit lower
triangular one.

## Other changes

* `nu` has its own argument in the inference functions. It is a prefix of
  `null`, so R's partial matching would otherwise have sent the number of
  cosine terms to the null value.
* `summary()` reports why intervals are unavailable instead of omitting them
  silently.
* Confidence regions keep the geometry of their method: a convex polygon for
  the increment hull, an ellipsoid for the other five. `aersn_projection()`,
  `aersn_slice()`, `aersn_contains()` and `plot()` work for both.
* New vignette, *Comparing inference methods on one estimate*.
* `sandwich` is a new dependency, used for the Andrews and Newey-West
  bandwidth rules.

## Scope

The structural-break test of Hong, Linton, McCabe, Sun and Wang (2024) and
the autocorrelation tests of Sun, Zhu and Linton (2025) remain outside the
package. No method here is an adjusted-range version of either: substituting
a different normalizer into a published statistic would need its own theory
and reference law.

# aersn 0.1.1

Review fixes. See `dev/FIXES_0.1.1.md` for the full list: reference-cache
keys, missing-value handling in `aersn_lm()`, degenerate profiles and failed
reference draws, quantile-table and parameter-selection validation, one-sided
tests, Monte Carlo uncertainty, numerical scaling, and plotting. Authorship
was completed to the five manuscript authors with Jiajing Sun as maintainer.

# aersn 0.1.0

Initial development version.

* Core construction: centered influence path (`aersn_path()`), increment-hull
  self-normalizer (`aersn_hull()`), gauge by linear programming
  (`aersn_gauge()`), support function (`aersn_support()`).
* Inference: `aersn()`, `aersn_test()`, `confint()` with simultaneous, joint
  and marginal intervals, `aersn_contrast()`, `aersn_region()`,
  `aersn_contains()`, two-dimensional projections and slices, plotting.
* Reference distributions: closed-form scalar law, matched-grid Monte Carlo
  references with reproducible seeds, batch Monte Carlo standard errors,
  caching, and the verified registry `aersn_registry`.
* Model interfaces: `aersn_mean()`, `aersn_lm()`, `aersn_gmm()`,
  `aersn_mle()`, and `aersn_target()`.
