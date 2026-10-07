#' Closed-form scalar adjusted-range reference law
#'
#' Distribution function, quantile function and density of the
#' continuous-path scalar reference law.  Let \eqn{Z \sim N(0, 1)} and let
#' \eqn{R} be the range (supremum minus infimum) of an independent standard
#' Brownian bridge on \eqn{[0, 1]}.  The signed ratio \eqn{M = Z / R} is the
#' limit of the signed scalar statistic, and \eqn{|M|} is the limit of the
#' two-sided statistic \eqn{T_n(\theta_0)} for \eqn{q = 1} (manuscript
#' Theorem 1 with \eqn{q = 1}, equation (scalar-reduction)).
#'
#' @details
#' The Supplement (Proposition "Density and distribution") gives
#' \deqn{f_M(u) = \tfrac12 - 12u^2\sum_{k \ge 1}\frac{k^2}{(u^2 + 4k^2)^{5/2}},
#'   \qquad
#'   \Pr(|M| \le c) = c - 2c^3\sum_{k \ge 1}(c^2 + 4k^2)^{-3/2},}
#' with truncation errors at most \eqn{c^3/(8N^2)} for the distribution
#' function and \eqn{3U^2/(16N^2)} for the density on \eqn{|u| \le U} after
#' \eqn{N} terms.  Poisson summation gives the exponentially convergent
#' tail representation
#' \deqn{\Pr(|M| > c) = 2\pi c^2\sum_{k \ge 1} k\,K_1(\pi k c),}
#' where \eqn{K_1} is the modified Bessel function of the second kind
#' (Supplement, equation (bessel-tail)), and differentiating it gives
#' \eqn{f_M(c) = \pi^2 c^2\sum_k k^2 K_0(\pi k c) - \pi c\sum_k k K_1(\pi k c)}.
#' The `"auto"` method uses the Bessel representations except very close
#' to zero, where the direct series is used; the two representations are
#' cross-checked in the package tests.
#'
#' The quantile function inverts the distribution function of \eqn{|M|}
#' with [stats::uniroot()].  Continuous-path critical values at two-sided
#' levels 10, 5 and 1 percent are 1.397390, 1.705776 and 2.367378
#' (Supplement Table "Continuous-path scalar critical values").
#'
#' These are continuous-path laws.  For a sample of size `n`, the
#' manuscript uses matched-grid quantiles, which are larger; see
#' [aersn_reference()] with `grid = "matched"`.
#'
#' @param c Nonnegative evaluation points for the distribution function of
#'   \eqn{|M|}.
#' @param u Evaluation points for the density of the signed ratio \eqn{M}.
#' @param p Probabilities strictly between 0 and 1.
#' @param method `"auto"`, `"bessel"` or `"series"`.
#' @param tol Truncation tolerance for the series representations.
#'
#' @return Numeric vectors: probabilities \eqn{\Pr(|M| \le c)}, quantiles of
#'   \eqn{|M|}, or density values of \eqn{M}.
#' @seealso [aersn_reference()].
#' @examples
#' aersn_scalar_cdf(c(1.397390, 1.705776, 2.367378))
#' aersn_scalar_quantile(c(0.90, 0.95, 0.99))
#' integrate(aersn_scalar_density, -Inf, Inf)$value   # integrates to one
#' @export
aersn_scalar_cdf <- function(c, method = c("auto", "bessel", "series"),
                             tol = 1e-13) {
  method <- match.arg(method)
  check_tolerance(tol)
  c <- as.numeric(c)
  if (any(!is.finite(c)) || any(c < 0)) {
    aersn_abort("`c` must contain finite nonnegative values.")
  }
  vapply(c, function(x) {
    if (x == 0) return(0)
    use_bessel <- switch(method, auto = x >= 0.02, bessel = TRUE,
                         series = FALSE)
    if (use_bessel) 1 - scalar_tail_bessel(x, tol) else
      scalar_cdf_series(x, tol)
  }, numeric(1L))
}

scalar_cdf_series <- function(x, tol) {
  N <- max(10L, ceiling(sqrt(x^3 / (8 * tol))))
  check_series_length(N, 1e7)
  k2 <- 4 * seq_len(N)^2
  x - 2 * x^3 * sum((x^2 + k2)^(-1.5))
}

scalar_tail_bessel <- function(x, tol) {
  if (x > 250) return(0) # below double-precision range
  ## Terms k K_1(pi k x) decay like exp(-pi k x); choose the number of
  ## terms from the target tolerance.
  K <- max(5L, ceiling((log(1 / tol) + 8) / (pi * x)))
  check_series_length(K, 2e6)
  k <- seq_len(K)
  2 * pi * x^2 * sum(k * besselK(pi * k * x, nu = 1))
}

## Avoid cancellation in 1 - F(x) for small upper-tail probabilities.
scalar_tail <- function(x, tol = 1e-13) {
  if (x == 0) return(1)
  if (x < 0.02) return(1 - scalar_cdf_series(x, tol))
  scalar_tail_bessel(x, tol)
}

check_tolerance <- function(tol) {
  check_scalar_number(tol, "tol", lower = 0, upper = 1)
  if (tol == 0 || tol == 1) aersn_abort("`tol` must lie strictly between 0 and 1.")
}

check_series_length <- function(n, limit) {
  if (!is.finite(n) || n > limit) {
    aersn_abort("The requested series needs too many terms for this argument ",
                "and tolerance. Use method = \"auto\" or relax `tol`.")
  }
}

#' @rdname aersn_scalar_cdf
#' @export
aersn_scalar_quantile <- function(p, tol = 1e-12) {
  check_tolerance(tol)
  p <- as.numeric(p)
  if (any(!is.finite(p)) || any(p <= 0) || any(p >= 1)) {
    aersn_abort("`p` must contain probabilities strictly between 0 and 1.")
  }
  vapply(p, function(pp) {
    upper <- 4
    while (aersn_scalar_cdf(upper) < pp) upper <- 2 * upper
    stats::uniroot(function(x) aersn_scalar_cdf(x) - pp,
                   interval = c(0, upper), tol = tol, maxiter = 10000L)$root
  }, numeric(1L))
}

#' @rdname aersn_scalar_cdf
#' @export
aersn_scalar_density <- function(u, method = c("auto", "bessel", "series"),
                                 tol = 1e-11) {
  method <- match.arg(method)
  check_tolerance(tol)
  u <- as.numeric(u)
  if (any(!is.finite(u))) aersn_abort("`u` must be finite.")
  vapply(abs(u), function(x) {
    use_bessel <- switch(method, auto = x >= 0.05, bessel = TRUE,
                         series = FALSE)
    if (use_bessel && x > 0) scalar_density_bessel(x, tol) else
      scalar_density_series(x, tol)
  }, numeric(1L))
}

scalar_density_series <- function(x, tol) {
  N <- max(10L, ceiling(sqrt(3 * max(x, 1e-3)^2 / (16 * tol))))
  check_series_length(N, 1e7)
  k <- seq_len(N)
  0.5 - 12 * x^2 * sum(k^2 / (x^2 + 4 * k^2)^2.5)
}

scalar_density_bessel <- function(x, tol) {
  if (x > 250) return(0)
  K <- max(5L, ceiling((log(1 / tol) + 8) / (pi * x)))
  check_series_length(K, 2e6)
  k <- seq_len(K)
  a <- pi * k * x
  pi^2 * x^2 * sum(k^2 * besselK(a, nu = 0)) -
    pi * x * sum(k * besselK(a, nu = 1))
}
