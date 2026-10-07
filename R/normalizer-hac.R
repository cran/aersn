## Kernel (HAC) long-run covariance estimation.
##
## HAC estimation is a covariance-estimation method, not a self-normalization
## method: the long-run covariance matrix is estimated consistently and the
## Wald statistic has a chi-squared limit. It is included so that the
## adjusted-range methods can be compared with the standard practice on the
## same data, tuning values and level.

## Kernel weights w(x), matching sandwich::kweights.
aersn_kernel_weight <- function(x, kernel) {
  switch(kernel,
         "Bartlett" = pmax(0, 1 - abs(x)),
         "Parzen" = {
           ax <- abs(x)
           ifelse(ax <= 0.5, 1 - 6 * ax^2 + 6 * ax^3,
                  ifelse(ax <= 1, 2 * (1 - ax)^3, 0))
         },
         "Quadratic Spectral" = {
           y <- 6 * pi * x / 5
           out <- numeric(length(y))
           small <- abs(y) < 1e-3
           y2 <- y[small]^2
           out[small] <- 1 - y2 / 10 + y2^2 / 280 - y2^3 / 15120
           yy <- y[!small]
           out[!small] <- 3 / yy^2 * (sin(yy) / yy - cos(yy))
           out
         },
         aersn_abort("Unsupported kernel '", kernel, "'."))
}

aersn_hac_kernels <- function() {
  c("Bartlett", "Parzen", "Quadratic Spectral")
}

aersn_hac_rules <- function() {
  c("short", "long", "andrews", "newey-west")
}

## Weighted sum of sample autocovariances.  `divisor` is the sample size used
## in the scaling: with VAR(1) prewhitening there are n - 1 residuals but the
## divisor stays n, which is the convention of sandwich::vcovHAC().
hac_matrix_from_weights <- function(E, w, divisor = nrow(E)) {
  m <- nrow(E)
  V <- crossprod(E) / divisor
  nz <- which(w != 0)
  for (l in nz) {
    G <- crossprod(E[seq_len(m - l), , drop = FALSE],
                   E[(l + 1L):m, , drop = FALSE]) / divisor
    V <- V + w[l] * (G + t(G))
  }
  (V + t(V)) / 2
}

## Andrews-Monahan VAR(1) prewhitening.
hac_prewhiten <- function(E) {
  n <- nrow(E)
  q <- ncol(E)
  units <- sqrt(colMeans(E^2))
  if (any(!is.finite(units)) || any(units <= 0)) {
    aersn_abort("VAR(1) prewhitening needs nonzero finite contribution scales.")
  }
  Es <- sweep(E, 2L, units, "/")
  X <- Es[-n, , drop = FALSE]
  Y <- Es[-1L, , drop = FALSE]
  fit <- qr(X)
  if (fit$rank < q) {
    aersn_abort("VAR(1) prewhitening failed: the lagged contribution matrix ",
                "is numerically singular. Use prewhite = FALSE.")
  }
  As <- t(qr.coef(fit, Y))
  resid <- sweep(Y - X %*% t(As), 2L, units, "*")
  Ds <- diag(q) - As
  evd <- svd(Ds, nu = 0L, nv = 0L)$d
  if (min(evd) <= 1e-10 * max(evd)) {
    aersn_abort("VAR(1) prewhitening failed: I - A is numerically singular, ",
                "so the recoloring step is not defined. Use prewhite = FALSE.")
  }
  A <- sweep(sweep(As, 1L, units, "*"), 2L, units, "/")
  D <- diag(q) - A
  Dinv <- sweep(sweep(solve(Ds), 1L, units, "*"), 2L, units, "/")
  list(resid = resid, A = A, D = D, Dinv = Dinv)
}

#' HAC long-run covariance estimate
#'
#' Estimates the long-run covariance matrix of the influence contributions of
#' an [aersn()] fit by kernel smoothing of the sample autocovariances, and
#' returns it as an [aersn_normalizer()] so that tests, intervals and regions
#' can be built from it exactly as for the other methods.
#'
#' @details
#' Write \eqn{e_t} for the influence contributions, demeaned when
#' `center = TRUE`.  With kernel weights \eqn{w_\ell = w(\ell/h)} and
#' bandwidth \eqn{h},
#' \deqn{\hat\Omega = \hat\Gamma_0 + \sum_{\ell = 1}^{n-1} w_\ell
#'   (\hat\Gamma_\ell + \hat\Gamma_\ell^\top), \qquad
#'   \hat\Gamma_\ell = \frac 1n \sum_{t=1}^{n-\ell} e_t e_{t+\ell}^\top .}
#' The statistic for \eqn{H_0: \theta_0 = v} is
#' \eqn{z^\top\hat\Omega^{-1}z} with \eqn{z = \sqrt n(\hat\theta_n - v)}, and
#' the reference law is \eqn{\chi^2_q} (or \eqn{\chi^2_m} for an
#' \eqn{m}-dimensional linear restriction).  All \eqn{n-1} lags are used; the
#' Bartlett and Parzen weights vanish beyond the bandwidth, while the
#' quadratic spectral kernel has unbounded support and is not truncated.
#'
#' **Kernels.** `"Bartlett"` gives \eqn{w(x) = 1 - |x|} for \eqn{|x| \le 1};
#' `"Parzen"` gives \eqn{1 - 6x^2 + 6|x|^3} for \eqn{|x| \le 1/2} and
#' \eqn{2(1 - |x|)^3} for \eqn{1/2 < |x| \le 1}; `"Quadratic Spectral"` gives
#' \eqn{3\{\sin(y)/y - \cos(y)\}/y^2} with \eqn{y = 6\pi x/5}.  These match
#' [sandwich::kweights()].
#'
#' **Bandwidth.** Distinguish three quantities that are easy to confuse.  The
#' real-valued bandwidth \eqn{h} enters through \eqn{w(\ell/h)}.  A lag
#' truncation \eqn{L} for the Bartlett kernel corresponds to \eqn{h = L + 1},
#' since \eqn{w_\ell = 1 - \ell/(L+1)} is then zero at \eqn{\ell = L+1}.  The
#' fixed-b fraction \eqn{b} of [aersn_fixed_b_normalizer()] is a different
#' parameterisation again.  `bandwidth` accepts:
#' * a positive number, used directly as \eqn{h};
#' * `"short"`: \eqn{L = \lfloor 4(n/100)^{2/9}\rfloor}, \eqn{h = L + 1}, the
#'   rule used for the Newey-West comparison in the manuscript;
#' * `"long"`: \eqn{L = \lceil 1.3\sqrt n\,\rceil}, \eqn{h = L + 1};
#' * `"andrews"`: the plug-in bandwidth of Andrews (1991) for the chosen
#'   kernel, from [sandwich::bwAndrews()] with equal weights on the
#'   coordinates.  For the Parzen and quadratic spectral kernels the returned
#'   bandwidths satisfy \eqn{h_{\mathrm{Parzen}} = h_{\mathrm{QS}}\times
#'   2.6614/1.3221} exactly, the ratio of Andrews' kernel constants; the
#'   manuscript's replication code uses that identity and the package checks
#'   it in its tests.
#' * `"newey-west"`: the lag-selection rule of Newey and West (1994), from
#'   [sandwich::bwNeweyWest()], applied as \eqn{L = \lfloor bw \rfloor} and
#'   \eqn{h = L + 1}.  This is the convention of [sandwich::NeweyWest()],
#'   which applies Bartlett weights; it is therefore supported for the
#'   Bartlett kernel only, and other kernels raise an error rather than
#'   silently using a different rule.
#'
#' Alternatively, `lag` sets the Bartlett or Parzen lag truncation directly,
#' giving \eqn{h = L + 1}; it is not available for the quadratic spectral
#' kernel, whose weights do not vanish at any finite lag.
#'
#' **Prewhitening.** `prewhite = TRUE` applies the Andrews-Monahan procedure:
#' a VAR(1) is fitted to the contributions by least squares, the kernel
#' estimate is formed from its residuals, and the result is recolored as
#' \eqn{(I - \hat A)^{-1}\hat\Omega_e (I - \hat A)^{-\top}}.  Bandwidth
#' selection then uses the residuals, and the scaling divisor stays \eqn{n}
#' although there are \eqn{n - 1} residuals, which is the convention of
#' [sandwich::vcovHAC()].  This is VAR prewhitening of the HAC
#' estimator and is unrelated to the lag-zero LDL partial prewhitening of
#' [aersn_ldl_normalizer()].  The default `FALSE` reproduces the manuscript's
#' comparison.
#'
#' **Finite-sample adjustment.** `adjust = TRUE` multiplies the estimate by
#' \eqn{n/(n-q)}, where \eqn{q} is the target dimension, not the total
#' number of fitted coefficients. For a selected target this can differ
#' from [sandwich::meatHAC()] applied to the full model.
#' The default `FALSE` reproduces the manuscript's
#' comparison.  The influence contributions already absorb the estimation of
#' nuisance parameters, so this multiplier is a finite-sample convention
#' rather than a correction derived for this framework.
#'
#' @param object An [aersn()] fit.
#' @param kernel `"Bartlett"`, `"Parzen"` or `"Quadratic Spectral"`.
#' @param bandwidth A positive number, or one of `"short"`, `"long"`,
#'   `"andrews"`, `"newey-west"`.  Ignored when `lag` is supplied.
#' @param lag Optional nonnegative integer lag truncation for the Bartlett or
#'   Parzen kernel, giving bandwidth `lag + 1`.
#' @param prewhite Apply VAR(1) prewhitening and recoloring.
#' @param adjust Apply the \eqn{n/(n-q)} finite-sample multiplier.
#' @param center Demean the influence contributions before forming the
#'   autocovariances.  For fitted contributions the column means are already
#'   zero and this has no effect.
#'
#' @return An [aersn_normalizer()] object of class
#'   `"aersn_normalizer_hac"`.  `tuning` records the kernel, the bandwidth
#'   rule, the realized numeric bandwidth, the implied Bartlett lag truncation
#'   where that is meaningful, and the prewhitening, adjustment and centering
#'   flags.
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. *Econometrica*, 59(3), 817-858.
#'
#' Newey, W. K. and West, K. D. (1994). Automatic lag selection in covariance
#' matrix estimation. *Review of Economic Studies*, 61(4), 631-653.
#'
#' Andrews, D. W. K. and Monahan, J. C. (1992). An improved
#' heteroskedasticity and autocorrelation consistent covariance matrix
#' estimator. *Econometrica*, 60(4), 953-966.
#' @seealso [aersn_normalizer()], [aersn_fixed_b_normalizer()],
#'   [aersn_ewc_lrv()], [aersn_compare()].
#' @examples
#' set.seed(4)
#' n <- 300
#' e <- matrix(rnorm(2 * n), n, 2)
#' Y <- e
#' for (t in 2:n) Y[t, ] <- 0.5 * Y[t - 1, ] + e[t, ]
#' fit <- aersn_mean(Y)
#' aersn_hac_lrv(fit, kernel = "Bartlett", bandwidth = "short")
#' aersn_hac_lrv(fit, kernel = "Quadratic Spectral", bandwidth = "andrews")
#' aersn_test(fit, method = "hac", kernel = "Parzen", bandwidth = "andrews")
#' @export
aersn_hac_lrv <- function(object, kernel = c("Bartlett", "Parzen",
                                             "Quadratic Spectral"),
                          bandwidth = "short", lag = NULL,
                          prewhite = FALSE, adjust = FALSE, center = TRUE) {
  check_aersn(object)
  kernel <- match.arg(kernel)
  check_flag(prewhite, "prewhite")
  check_flag(adjust, "adjust")
  check_flag(center, "center")
  n <- object$n
  q <- object$q
  E <- object$psi
  if (center) E <- sweep(E, 2L, colMeans(E), "-")
  A <- NULL
  if (prewhite) {
    if (n < q + 2L) {
      aersn_abort("VAR(1) prewhitening needs more than q + 1 observations.")
    }
    pw <- hac_prewhiten(E)
    E_used <- pw$resid
    A <- pw$A
  } else {
    E_used <- E
  }
  bw <- hac_bandwidth(E_used, kernel, bandwidth, lag, n = nrow(E_used), q = q)
  m <- nrow(E_used)
  w <- aersn_kernel_weight(seq_len(m - 1L) / bw$h, kernel)
  V <- hac_matrix_from_weights(E_used, w, divisor = n)
  if (prewhite) {
    Dinv <- pw$Dinv
    V <- Dinv %*% V %*% t(Dinv)
    V <- (V + t(V)) / 2
  }
  if (adjust) {
    if (n <= q) aersn_abort("`adjust = TRUE` needs n > q.")
    V <- V * n / (n - q)
  }
  dimnames(V) <- list(object$names, object$names)
  new_normalizer(
    "hac", q = q, n = n, names = object$names, matrix = V,
    tuning = list(kernel = kernel, bandwidth_rule = bw$rule,
                  bandwidth = bw$h,
                  lag_truncation = bw$lag,
                  prewhite = prewhite, adjust = adjust, center = center),
    diagnostics = list(
      matrix_condition = matrix_condition(V),
      nonzero_lags = sum(w != 0),
      max_abs_weight_at_last_lag = abs(w[length(w)]),
      var1_max_abs_coefficient = if (is.null(A)) NA_real_ else max(abs(A))),
    reference_args = list(df = q),
    centering = object$centering, nodes = object$nodes, source = object,
    notes = "Consistent covariance estimation with a chi-squared reference; not a self-normalization method.")
}

## Resolve the bandwidth, returning the numeric h, the rule label, and the
## implied Bartlett/Parzen lag truncation when one is meaningful.
hac_bandwidth <- function(E, kernel, bandwidth, lag, n, q) {
  if (!is.null(lag)) {
    if (kernel == "Quadratic Spectral") {
      aersn_abort("`lag` sets a lag truncation, which is not defined for the ",
                  "quadratic spectral kernel: its weights do not vanish at ",
                  "any finite lag. Supply a numeric `bandwidth` instead.")
    }
    check_scalar_number(lag, "lag", lower = 0, integer = TRUE)
    return(list(h = as.numeric(lag) + 1, rule = "user lag truncation",
                lag = as.integer(lag)))
  }
  if (is.numeric(bandwidth)) {
    check_scalar_number(bandwidth, "bandwidth", lower = 0)
    if (bandwidth <= 0) aersn_abort("`bandwidth` must be positive.")
    return(list(h = as.numeric(bandwidth), rule = "user bandwidth",
                lag = if (kernel == "Quadratic Spectral") NA_integer_ else
                  as.integer(max(0, ceiling(bandwidth) - 1L))))
  }
  if (!is.character(bandwidth) || length(bandwidth) != 1L) {
    aersn_abort("`bandwidth` must be a positive number or one of ",
                paste0("\"", aersn_hac_rules(), "\"", collapse = ", "), ".")
  }
  rule <- match.arg(bandwidth, aersn_hac_rules())
  if (rule == "short") {
    L <- floor(4 * (n / 100)^(2 / 9))
    return(list(h = L + 1, rule = "short (floor(4 (n/100)^(2/9)))",
                lag = as.integer(L)))
  }
  if (rule == "long") {
    L <- ceiling(1.3 * sqrt(n))
    return(list(h = L + 1, rule = "long (ceiling(1.3 sqrt(n)))",
                lag = as.integer(L)))
  }
  if (rule == "andrews") {
    h <- sandwich::bwAndrews(E, kernel = kernel, weights = rep(1, q),
                             prewhite = FALSE)
    if (!is.finite(h) || h <= 0) {
      aersn_abort("The Andrews bandwidth rule returned ", h,
                  "; supply a numeric bandwidth instead.")
    }
    return(list(h = as.numeric(h), rule = "Andrews (1991) plug-in",
                lag = if (kernel == "Quadratic Spectral") NA_integer_ else
                  as.integer(max(0, ceiling(h) - 1L))))
  }
  ## newey-west
  if (kernel != "Bartlett") {
    aersn_abort("`bandwidth = \"newey-west\"` follows the convention of ",
                "sandwich::NeweyWest(), which applies Bartlett weights to ",
                "bandwidth floor(bw) + 1. It is supported for the ",
                "Bartlett kernel only. For the ", kernel, " kernel use ",
                "bandwidth = \"andrews\", a numeric bandwidth, or, for ",
                "Parzen, an explicit lag.")
  }
  L <- floor(sandwich::bwNeweyWest(E, kernel = kernel, weights = rep(1, q),
                                   prewhite = FALSE))
  if (!is.finite(L) || L < 0) {
    aersn_abort("The Newey-West lag rule returned ", L,
                "; supply a numeric bandwidth instead.")
  }
  list(h = L + 1, rule = "Newey-West (1994) automatic lag",
       lag = as.integer(L))
}
