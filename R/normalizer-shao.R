## Quadratic self-normalization (Shao 2010) from the centered influence path.

#' Quadratic self-normalizer
#'
#' Builds the quadratic self-normalizer of Shao (2010) from the centered
#' influence path of an [aersn()] fit.  The normalizer is the discretized
#' integral of the outer product of the path, and the statistic is the Wald
#' form with that matrix.  Like the increment hull it needs no bandwidth,
#' kernel or block length, and its reference law is nonstandard.
#'
#' @details
#' With the centered path \eqn{\hat G_n(k/n)}, \eqn{k = 0, \ldots, n}, the
#' calendar-time normalizer is
#' \deqn{V^{\mathrm Q}_n = \frac 1n \sum_{k=1}^{n}
#'   \hat G_n(k/n)\,\hat G_n(k/n)^\top ,}
#' the origin contributing zero.  The statistic for \eqn{H_0: \theta_0 = v} is
#' \eqn{z^\top (V^{\mathrm Q}_n)^{-1} z} with
#' \eqn{z = \sqrt n(\hat\theta_n - v)}, whose limit is
#' \eqn{B_q(1)^\top \{\int_0^1 \mathbb B_q \mathbb B_q^\top\}^{-1} B_q(1)}.
#' For \eqn{q = 1} the two-sided statistic is the square of
#' \eqn{|z| / \{n^{-1}\sum_k \hat G_n(k/n)^2\}^{1/2}}.
#'
#' **Integration rule.** Under nonuniform variance accumulation the
#' manuscript integrates with respect to the variance-accumulation profile
#' rather than calendar time.  With `integration = "profile"` the normalizer is
#' \deqn{V^{\mathrm Q,\tau}_n = \sum_{k=1}^{n} \hat G^\tau_n(k/n)\,
#'   \hat G^\tau_n(k/n)^\top \{\hat\tau_n(k/n) - \hat\tau_n((k-1)/n)\},}
#' the Riemann-Stieltjes sum with the estimated profile increments as weights.
#' With a linear profile the two rules coincide exactly.  The reference law is
#' simulated with the same rule and the same nodes, so a calendar-time
#' statistic is never combined with a profile-weighted reference.  The
#' adjusted-range methods do not need this correction, because projected
#' ranges are invariant under a monotone change of time.
#'
#' @param object An [aersn()] fit.
#' @param integration `"calendar"` for equal weights \eqn{1/n}, `"profile"`
#'   for the increments of the fit's variance-accumulation profile.  The
#'   default follows the fit: calendar-time centering gives `"calendar"`,
#'   profile centering gives `"profile"`.
#'
#' @return An [aersn_normalizer()] object of class
#'   `"aersn_normalizer_shao"`, carrying the matrix \eqn{V^{\mathrm Q}_n} and
#'   the integration weights actually used.
#' @references
#' Shao, X. (2010). A self-normalized approach to confidence interval
#' construction in time series. *Journal of the Royal Statistical Society:
#' Series B*, 72(3), 343-366. \doi{10.1111/j.1467-9868.2009.00737.x}
#' @seealso [aersn_normalizer()], [aersn_compare()].
#' @examples
#' set.seed(3)
#' fit <- aersn_mean(matrix(rnorm(400), 200, 2))
#' nz <- aersn_shao_normalizer(fit)
#' nz$matrix
#' aersn_test(fit, method = "shao", draws = 400, seed = 1)
#' @export
aersn_shao_normalizer <- function(object,
                                  integration = c("auto", "calendar", "profile")) {
  check_aersn(object)
  integration <- match.arg(integration)
  if (integration == "auto") {
    integration <- if (object$centering == "profile") "profile" else "calendar"
  }
  if (integration == "profile" && object$centering != "profile") {
    aersn_abort("`integration = \"profile\"` needs a fit built with a ",
                "variance-accumulation profile; this fit uses calendar-time ",
                "centering. See aersn_profile().")
  }
  n <- object$n
  q <- object$q
  G <- object$path$G[-1L, , drop = FALSE]     # k = 1, ..., n; the origin is zero
  w <- if (integration == "profile") diff(object$nodes) else rep(1 / n, n)
  if (abs(sum(w) - 1) > 1e-10) {
    aersn_abort("Internal: the integration weights do not sum to one.")
  }
  V <- crossprod(G * sqrt(w))
  dimnames(V) <- list(object$names, object$names)
  new_normalizer(
    "shao", q = q, n = n, names = object$names, matrix = V,
    tuning = list(integration = integration),
    diagnostics = list(
      weight_min = min(w), weight_max = max(w),
      matrix_condition = matrix_condition(V)),
    reference_args = list(integration = integration),
    centering = object$centering, nodes = object$nodes, source = object,
    notes = if (integration == "profile")
      "Integration uses the estimated variance-accumulation profile increments." else
        character())
}

matrix_condition <- function(V) {
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) return(Inf)
  max(ev) / min(ev)
}
