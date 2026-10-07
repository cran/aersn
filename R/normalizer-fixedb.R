## Bartlett fixed-b long-run variance estimation.

## Bartlett fixed-b matrix from the centered path.  G has n + 1 rows, the
## first being the origin.
fixed_b_matrix <- function(G, b) {
  N <- nrow(G) - 1L
  fixed_b_matrix_from_m(G, as.integer(round(b * N)))
}

## The same estimator addressed by its integer lag truncation, so that a
## reference law is simulated with exactly the truncation used on the sample.
fixed_b_matrix_from_m <- function(G, m) {
  N <- nrow(G) - 1L
  m <- as.integer(m)
  if (m < 1L || m > N) {
    aersn_abort("The fixed-b lag truncation must lie in 1..", N, "; got ", m, ".")
  }
  b_grid <- m / N
  level <- crossprod(G[-1L, , drop = FALSE]) / N
  left <- G[seq_len(N - m + 1L), , drop = FALSE]
  right <- G[(m + 1L):(N + 1L), , drop = FALSE]
  lagged <- crossprod(left, right) / N
  V <- (2 / b_grid) * level - (lagged + t(lagged)) / b_grid
  list(V = (V + t(V)) / 2, m = m, b_grid = b_grid)
}

#' Bartlett fixed-b normalizer
#'
#' Estimates the long-run covariance matrix with the Bartlett kernel and a
#' bandwidth proportional to the sample size, \eqn{M = bn}, and returns it as
#' an [aersn_normalizer()].  Because the bandwidth fraction \eqn{b} is held
#' fixed as the sample grows, the estimate does not converge to the long-run
#' covariance matrix; the randomness that remains in the limit is carried by a
#' nonstandard reference law that depends on \eqn{b}.
#'
#' @details
#' The estimate is computed from the centered influence path through the
#' partial-sum identity
#' \deqn{\hat\Omega_b = \frac{2}{b_{\mathrm{grid}}}\cdot
#'   \frac 1n\sum_{k=1}^{n} \hat G_{n,k}\hat G_{n,k}^\top
#'   - \frac{1}{b_{\mathrm{grid}}}\cdot\frac 1n\sum_{k=0}^{n-m}
#'   \bigl(\hat G_{n,k}\hat G_{n,k+m}^\top
#'         + \hat G_{n,k+m}\hat G_{n,k}^\top\bigr),}
#' which is algebraically the Bartlett estimator with bandwidth
#' \eqn{m = \mathrm{round}(bn)} and weights \eqn{1-\ell/m}
#' for \eqn{\ell < m}, applied to the demeaned contributions.
#' In [aersn_hac_lrv()] this corresponds to `lag = m - 1`.
#'
#' **Requested and realized b.** The bandwidth in grid intervals must be an integer, so the
#' realized fraction is \eqn{b_{\mathrm{grid}} = m/n}, which differs from the
#' requested \eqn{b} by at most \eqn{1/(2n)}.  The realized value is what
#' enters the estimate, is recorded in `tuning$b_grid`, and is the value used
#' to key and simulate the reference law.  A request with
#' \eqn{\mathrm{round}(bn) = 0} is rejected rather than silently replaced.
#'
#' **Relation to quadratic self-normalization.** At \eqn{b = 1} we have
#' \eqn{m = n}, the lagged term involves only the two endpoints of the path,
#' which are zero, and the identity reduces to
#' \eqn{\hat\Omega_1 = 2V^{\mathrm Q}_n} exactly in finite samples, where
#' \eqn{V^{\mathrm Q}_n} is the [aersn_shao_normalizer()] matrix.  The Wald
#' statistic is therefore exactly half the quadratic self-normalized
#' statistic, and the fixed-b reference law at \eqn{b = 1} is exactly half the
#' quadratic reference law, so the two tests reject on the same samples.  The
#' package tests verify both the matrix identity and the equality of the
#' resulting decisions using common reference draws.
#'
#' **Other kernels.** Only the Bartlett kernel is supported here, because the
#' package supplies a matching fixed-b reference law only for that kernel.
#' A Parzen or quadratic spectral estimate with a fixed bandwidth fraction
#' would need its own reference law; use [aersn_hac_lrv()] with a numeric
#' bandwidth if a consistent chi-squared reference is intended instead.
#'
#' @param object An [aersn()] fit.
#' @param b The bandwidth fraction, in \eqn{(0, 1]}.
#'
#' @return An [aersn_normalizer()] object of class
#'   `"aersn_normalizer_fixedb"`, whose `tuning` records the requested `b`,
#'   the realized `b_grid`, and the integer bandwidth `m`.
#' @references
#' Kiefer, N. M. and Vogelsang, T. J. (2005). A new asymptotic theory for
#' heteroskedasticity-autocorrelation robust tests. *Econometric Theory*,
#' 21(6), 1130-1164.
#' @seealso [aersn_normalizer()], [aersn_shao_normalizer()],
#'   [aersn_hac_lrv()], [aersn_compare()].
#' @examples
#' set.seed(5)
#' fit <- aersn_mean(matrix(rnorm(400), 200, 2))
#' nz <- aersn_fixed_b_normalizer(fit, b = 0.5)
#' nz$tuning$b_grid
#' ## b = 1 gives exactly twice the quadratic self-normalizer
#' max(abs(aersn_fixed_b_normalizer(fit, b = 1)$matrix -
#'           2 * aersn_shao_normalizer(fit)$matrix))
#' @export
aersn_fixed_b_normalizer <- function(object, b = 0.5) {
  check_aersn(object)
  check_scalar_number(b, "b", lower = 0, upper = 1)
  if (b <= 0) aersn_abort("`b` must be strictly positive.")
  n <- object$n
  m <- as.integer(round(b * n))
  if (m < 1L) {
    aersn_abort("`b = ", b, "` gives a lag truncation round(b * n) = 0 at ",
                "n = ", n, ", for which the fixed-b estimator is not defined. ",
                "Choose b strictly above ",
                signif(0.5 / n, 3), " so that round(b * n) is at least one.")
  }
  fb <- fixed_b_matrix(object$path$G, b)
  V <- fb$V
  dimnames(V) <- list(object$names, object$names)
  new_normalizer(
    "fixedb", q = object$q, n = n, names = object$names, matrix = V,
    tuning = list(b = b, b_grid = fb$b_grid, m = fb$m, kernel = "Bartlett"),
    diagnostics = list(matrix_condition = matrix_condition(V),
                       b_rounding = fb$b_grid - b),
    reference_args = list(b_grid = fb$b_grid, m = fb$m),
    centering = object$centering, nodes = object$nodes, source = object,
    notes = paste0("Fixed-smoothing reference law for b = ",
                   fmt_num(fb$b_grid, 6),
                   "; the estimate is not consistent for the long-run covariance."))
}
