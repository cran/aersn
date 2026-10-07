#' aersn: Affine-Equivariant Adjusted-Range Self-Normalization
#'
#' The package implements affine-equivariant adjusted-range
#' self-normalization for inference on a fixed-dimensional parameter of a
#' dependent time series. The construction takes the centered partial-sum
#' path of estimated influence contributions, forms the convex hull of all
#' its increments (the *increment hull*), and normalizes the scaled
#' estimation error by the gauge of that hull. The statistic has a pivotal
#' limit that depends only on the dimension of the parameter of interest,
#' so no long-run covariance matrix, bandwidth, kernel or block length has
#' to be chosen.
#'
#' @section What the package computes:
#' Let \eqn{\hat\theta_n \in R^q} estimate \eqn{\theta_0} and let
#' \eqn{\hat\psi_{n,t} \in R^q}, \eqn{t = 1, \ldots, n}, be estimated
#' influence contributions with
#' \eqn{\sqrt n(\hat\theta_n - \theta_0) = n^{-1/2}\sum_t \psi_t + o_p(1)}.
#' The centered path on the observation grid is
#' \deqn{\hat G_n(k/n) = n^{-1/2}\{\sum_{t \le k}\hat\psi_{n,t}
#'   - \tau(k/n)\sum_{t \le n}\hat\psi_{n,t}\}, \quad k = 0, \ldots, n,}
#' with \eqn{\tau(r) = r} under calendar-time centering.  Its increment hull
#' is \eqn{K(\hat G_n) = \mathrm{conv}\{\hat G_n(r) - \hat G_n(s)\}}, whose
#' support function in direction \eqn{u} equals the adjusted range (maximum
#' minus minimum) of the projected path \eqn{u^\top \hat G_n}.  The test
#' statistic for \eqn{H_0: \theta_0 = v} is the gauge
#' \eqn{T_n(v) = \gamma_{K(\hat G_n)}\{\sqrt n(\hat\theta_n - v)\}},
#' evaluated by a linear program over the \eqn{n + 1} path values.  The
#' confidence region is \eqn{\hat\theta_n + (c/\sqrt n) K(\hat G_n)}, where
#' \eqn{c} is a quantile of the reference law
#' \eqn{\gamma_{K(B_q)}(Z_q)} formed from a standard Brownian bridge
#' \eqn{B_q} and an independent standard normal vector \eqn{Z_q}.
#'
#' @section Comparator methods:
#' Five further methods are available on the same fit, the same influence
#' contributions and the same level, through the `method` argument of
#' [aersn_test()], [confint.aersn()], [aersn_contrast()] and
#' [aersn_region()], and side by side through [aersn_compare()]:
#' componentwise adjusted ranges after lag-zero LDL partial prewhitening
#' ([aersn_ldl_normalizer()]), quadratic self-normalization
#' ([aersn_shao_normalizer()]), kernel long-run covariance estimation
#' ([aersn_hac_lrv()]), Bartlett fixed-b ([aersn_fixed_b_normalizer()]) and
#' the equal-weighted cosine method ([aersn_ewc_lrv()]).  Their statistics are
#' on different scales and are never compared as numbers; the comparable
#' quantities are p-values, decisions and interval widths.
#'
#' @section Main functions:
#' * [aersn()] builds the inference object from an estimate and its
#'   influence contributions; [aersn_mean()], [aersn_lm()], [aersn_gmm()]
#'   and [aersn_mle()] compute the contributions for common estimators and
#'   call the same core construction.
#' * [aersn_test()], [confint.aersn()], [aersn_contrast()],
#'   [aersn_region()] and [aersn_contains()] provide tests, simultaneous
#'   intervals, and confidence-region membership.
#' * [aersn_reference()] simulates or evaluates reference distributions;
#'   [aersn_scalar_cdf()] and [aersn_scalar_quantile()] give the closed-form
#'   scalar law.
#' * [aersn_path()], [aersn_hull()], [aersn_gauge()] and [aersn_support()]
#'   expose the geometric core, and [aersn_normalizer()] the normalizer of
#'   any method.
#' * [aersn_compare()] applies several methods to one estimate and reports
#'   the tuning values each of them used.
#'
#' @section Validity conditions:
#' Inference is asymptotically valid under the conditions of the
#' manuscript: (i) the partial-sum process of the influence contributions
#' satisfies a functional central limit theorem with a nonsingular scale
#' matrix; (ii) the estimator is asymptotically linear in those
#' contributions; and (iii) the feasible increment hull converges in
#' Hausdorff distance to the hull of the infeasible centered path.
#' The package cannot verify these conditions for arbitrary estimators;
#' see the vignette *Supplied influence contributions* and the
#' documentation of each model interface for what must hold.
#'
#' @references
#' Hong, Y., Lin, Z., Linton, O., Newey, W. K. and Sun, J. (2026).
#' Affine-equivariant adjusted-range self-normalization.
#' Cambridge Working Papers in Economics No. 2678;
#' Janeway Institute Working Paper No. 2637.
#' [Publication page](https://www.janeway.econ.cam.ac.uk/publication/affine-equivariant-adjusted-range-self-normalization).
#'
#' Hong, Y., Linton, O., McCabe, B., Sun, J. and Wang, S. (2024).
#' Kolmogorov-Smirnov type testing for structural breaks: a new
#' adjusted-range based self-normalization approach.
#' *Journal of Econometrics*, 238(2), 105603.
#' \doi{10.1016/j.jeconom.2023.105603}
#'
#' Shao, X. (2010). A self-normalized approach to confidence interval
#' construction in time series. *Journal of the Royal Statistical
#' Society: Series B*, 72(3), 343-366.
#' \doi{10.1111/j.1467-9868.2009.00737.x}
#'
#' @seealso [Working paper and publication details](https://www.janeway.econ.cam.ac.uk/publication/affine-equivariant-adjusted-range-self-normalization).
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats quantile confint rnorm uniroot pnorm qnorm sd
#' @importFrom stats model.matrix residuals coef lm setNames nobs na.action formula
#' @importFrom graphics plot polygon points abline segments par matplot
#' @importFrom graphics legend text title lines
#' @importFrom grDevices chull adjustcolor
#' @importFrom utils data head
## usethis namespace: end
NULL
