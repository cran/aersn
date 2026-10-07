## Equal-weighted cosine (EWC) heteroskedasticity- and autocorrelation-robust
## covariance estimation.

## Type II discrete cosine basis, orthonormal on the observation grid and
## orthogonal to the constant.
ewc_basis <- function(n, nu) {
  tt <- seq_len(n) - 0.5
  outer(tt, seq_len(nu), function(t, j) sqrt(2 / n) * cos(pi * j * t / n))
}

#' Equal-weighted cosine (EWC) covariance estimate
#'
#' Estimates the long-run covariance matrix by averaging the outer products of
#' the projections of the influence contributions on a cosine basis, and
#' returns it as an [aersn_normalizer()].  EWC is a heteroskedasticity- and
#' autocorrelation-robust (HAR) method in the fixed-smoothing class: the
#' reference law holds the number of basis terms fixed as the sample grows,
#' retaining the randomness of the covariance estimate through a scaled F
#' distribution. The default number of terms increases with sample size;
#' at each sample size the corresponding scaled F reference is used.
#'
#' @details
#' With \eqn{\nu} basis terms the type II cosine functions are
#' \deqn{\phi_j(t) = \sqrt{2/n}\,\cos\{\pi j (t - 1/2)/n\},
#'   \qquad j = 1, \ldots, \nu,\ t = 1, \ldots, n,}
#' which are orthonormal on the observation grid and exactly orthogonal to the
#' constant, so demeaning the contributions does not change the estimate.  The
#' projections and the estimate are
#' \deqn{\Lambda_j = \sum_{t=1}^n \phi_j(t)\, e_t, \qquad
#'   \hat\Omega_{\mathrm{EWC}} = \frac 1\nu \sum_{j=1}^\nu
#'     \Lambda_j \Lambda_j^\top .}
#' The statistic for \eqn{H_0: \theta_0 = v} is
#' \eqn{z^\top\hat\Omega_{\mathrm{EWC}}^{-1}z} with
#' \eqn{z = \sqrt n(\hat\theta_n - v)}.
#'
#' **Reference law.** For a target of dimension \eqn{m} the statistic is
#' referred to a scaled F law: \eqn{\{(\nu - m + 1)/(\nu m)\}\,T} is compared
#' with \eqn{F_{m,\,\nu - m + 1}}, equivalently the critical value is
#' \eqn{\nu m/(\nu - m + 1)} times the \eqn{F_{m,\nu-m+1}} quantile.  This is
#' the Hotelling form, and it requires \eqn{\nu \ge m}.  For a linear
#' restriction the dimension \eqn{m} of the restriction is used, not the full
#' parameter dimension.  The law is the fixed-smoothing asymptotic reference
#' for dependent data; it is not an exact finite-sample distribution.
#'
#' **Number of terms.** The default is
#' \eqn{\nu = \lfloor 0.4\, n^{2/3}\rfloor}, the rule used in the
#' manuscript's comparison, following the recommendation of
#' Lazarus, Lewis, Stock and Watson (2018).  Admissible values satisfy
#' \eqn{q \le \nu \le n - 1}; the upper limit is the number of cosine terms
#' that remain orthonormal on the grid.  Larger \eqn{\nu} reduces the variance
#' of the estimate and moves the reference law towards chi-squared, at the
#' cost of a larger bias under strong dependence.
#'
#' @param object An [aersn()] fit.
#' @param nu Number of cosine terms.  The default applies the rule above.
#' @param center Demean the influence contributions first.  The estimate does
#'   not depend on this because the basis is orthogonal to the constant; the
#'   argument exists so that the convention can be stated explicitly.
#'
#' @return An [aersn_normalizer()] object of class
#'   `"aersn_normalizer_ewc"`, whose `tuning` records `nu`, the rule used and
#'   the centering flag.
#' @references
#' Lazarus, E., Lewis, D. J., Stock, J. H. and Watson, M. W. (2018).
#' HAR inference: recommendations for practice. *Journal of Business and
#' Economic Statistics*, 36(4), 541-559.
#' @seealso [aersn_normalizer()], [aersn_hac_lrv()],
#'   [aersn_fixed_b_normalizer()], [aersn_compare()].
#' @examples
#' set.seed(6)
#' fit <- aersn_mean(matrix(rnorm(600), 300, 2))
#' nz <- aersn_ewc_lrv(fit)
#' nz$tuning$nu
#' aersn_test(fit, method = "ewc")
#' aersn_test(fit, method = "ewc", nu = 20)
#' @export
aersn_ewc_lrv <- function(object, nu = NULL, center = TRUE) {
  check_aersn(object)
  check_flag(center, "center")
  n <- object$n
  q <- object$q
  rule <- "floor(0.4 n^(2/3))"
  if (is.null(nu)) {
    nu <- floor(0.4 * n^(2 / 3))
  } else {
    check_scalar_number(nu, "nu", lower = 1, integer = TRUE)
    rule <- "user supplied"
  }
  nu <- as.integer(nu)
  if (nu < q) {
    aersn_abort("EWC needs at least as many cosine terms as parameters: ",
                "nu = ", nu, " but q = ", q, ". The reference law has ",
                "nu - q + 1 denominator degrees of freedom, and the estimate ",
                "is singular when nu < q.")
  }
  if (nu > n - 1L) {
    aersn_abort("EWC accepts at most n - 1 = ", n - 1L, " cosine terms on a ",
                "grid of n = ", n, " observations; nu = ", nu,
                " would repeat or annihilate basis vectors.")
  }
  E <- object$psi
  if (center) E <- sweep(E, 2L, colMeans(E), "-")
  Phi <- ewc_basis(n, nu)
  Lambda <- crossprod(Phi, E)          # nu x q
  V <- crossprod(Lambda) / nu
  V <- (V + t(V)) / 2
  dimnames(V) <- list(object$names, object$names)
  new_normalizer(
    "ewc", q = q, n = n, names = object$names, matrix = V,
    tuning = list(nu = nu, rule = rule, center = center,
                  basis = "type II cosine"),
    diagnostics = list(matrix_condition = matrix_condition(V),
                       basis_max_abs_column_sum = max(abs(colSums(Phi))),
                       smoothing_fraction = nu / n),
    reference_args = list(nu = nu),
    centering = object$centering, nodes = object$nodes, source = object,
    notes = "Scaled F reference at the selected number of cosine terms; its fixed-smoothing derivation holds nu fixed.")
}
