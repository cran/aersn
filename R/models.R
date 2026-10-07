## Model interfaces.  Each interface computes the estimate and the
## observation-level influence contributions of the manuscript and then
## calls aersn(); none of them re-implements the geometric construction.

#' Sample mean of a scalar or vector time series
#'
#' Adjusted-range inference for the mean.  The estimate is the sample
#' mean and the influence contributions are the demeaned observations,
#' \eqn{\hat\psi_{n,t} = Y_t - \bar Y_n} (manuscript Section 2.1).  The
#' centered path is the scaled cumulative sum of the demeaned series.
#'
#' @param Y Numeric vector (scalar series) or `n x q` matrix with
#'   observations in rows, in time order.
#' @param profile Optional variance-accumulation profile; see
#'   [aersn_profile()].  Because the demeaned observations sum to zero,
#'   a profile does not change the path; it changes the reference grid.
#'   The manuscript shows that this is not the required profile-centered
#'   path for the sample mean unless the profile is linear, so a profile
#'   should be supplied here only when its use is justified separately.
#' @param names Optional parameter names.
#' @return An [aersn()] object.
#' @details Validity requires the functional central limit theorem for the
#'   partial sums of \eqn{Y_t - \theta_0} with nonsingular long-run
#'   covariance (manuscript Assumption 1(i)); the representation is exact
#'   and the feasible path equals the infeasible centered path, so
#'   Assumption 1(ii)-(iii) hold automatically.
#' @examples
#' set.seed(11)
#' x <- arima.sim(list(ar = 0.5), 300)
#' fit <- aersn_mean(x)
#' aersn_test(fit, 0, reference = "continuous")
#' confint(fit, draws = 5000)
#' @export
aersn_mean <- function(Y, profile = NULL, names = NULL) {
  Y <- as_contribution_matrix(Y, "Y")
  est <- colMeans(Y)
  psi <- sweep(Y, 2L, est)
  names <- parameter_names(names, ncol(Y), psi = Y)
  aersn(est, psi, profile = profile, names = names,
        model = list(type = "mean",
                     description = "sample mean (psi_t = Y_t - Ybar)"),
        call = match.call())
}

#' Smooth generalized method of moments
#'
#' Adjusted-range inference for a target \eqn{\theta = h(\beta)} of a
#' GMM estimator \eqn{\hat\beta_n} with moment contributions
#' \eqn{m_t(\beta) \in R^{d_m}}, \eqn{d_m \ge p}.  With
#' \eqn{\hat D = n^{-1}\sum_t \partial m_t(\hat\beta_n)/\partial\beta^\top},
#' weighting matrix \eqn{\hat W} and
#' \eqn{\hat M = (\hat D^\top \hat W \hat D)^{-1}\hat D^\top \hat W}, the
#' influence contributions are
#' \deqn{\hat\psi_{n,t} = -\hat{\dot h}_n \hat M m_t(\hat\beta_n),}
#' manuscript Proposition "Feasible influence paths for smooth GMM".
#' The matrix \eqn{\hat M} uses the Jacobian with respect to all
#' coordinates of \eqn{\beta}, so the joint estimation of nuisance
#' coefficients is retained in the contributions of the target.
#'
#' @details
#' The estimator must (approximately) solve the first-order condition
#' \eqn{\hat D^\top \hat W \bar m_n(\hat\beta_n) = o_p(n^{-1/2})}; the
#' function reports each component of this scaled condition divided by
#' the sample standard deviation of the corresponding component of
#' \eqn{\hat D^\top\hat W m_t}. It warns when the largest absolute
#' standardized residual exceeds `.Machine$double.eps^0.25` (about
#' `1.2e-4`). This numerical convergence check is not a statistical test
#' or a verification of the asymptotic assumptions. The
#' sufficient conditions are manuscript Assumption "Sufficient conditions
#' for smooth GMM": stationarity and ergodicity of the moment process,
#' \eqn{\sqrt n}-consistency, full column rank of \eqn{D_m}, consistent
#' Jacobian and weighting matrix, smoothness of the moments, and a joint
#' functional central limit theorem for the moments and their Jacobian.
#' Under these conditions the feasible path is uniformly equivalent to the
#' path based on the true influence function.  This includes OLS, IV, and
#' stacked local-projection systems; it does not cover non-smooth moments
#' (for example quantile regression) or estimators with nonstandard rates.
#'
#' @param moments `n x d_m` matrix of moment contributions
#'   \eqn{m_t(\hat\beta_n)} evaluated at the estimate, in time order.
#' @param jacobian `d_m x p` matrix \eqn{\hat D}, the average derivative of
#'   the moments with respect to \eqn{\beta} at the estimate.
#' @param estimate The estimate \eqn{\hat\beta_n} (length `p`).
#' @param weight `d_m x d_m` positive-definite weighting matrix; `NULL`
#'   gives the identity (irrelevant when `d_m = p`).
#' @param target The parameters of interest: `NULL` for all of `beta`;
#'   integer indices or names selecting coordinates; a `q x p` matrix of
#'   full row rank (linear targets); or a list with elements `h` (a function
#'   of `beta`) and `jacobian` (the `q x p` derivative, a matrix or a
#'   function of `beta`).
#' @param profile Optional variance-accumulation profile; see
#'   [aersn_profile()].  For GMM the manuscript verifies feasibility only
#'   under calendar-time centering (\eqn{\tau(r) = r}).
#' @param names Optional names for the target parameters.
#' @return An [aersn()] object for the target, with `model$foc` holding the
#'   scaled first-order condition \eqn{\sqrt n \hat D^\top\hat W\bar m_n}.
#'   `model$foc_standardized` and `model$foc_tolerance` record the
#'   standardized residuals and the numerical warning threshold.
#' @examples
#' ## Instrumental-variables example (synthetic): y = 1 + 0.5 x + u,
#' ## x endogenous, instrument z.
#' set.seed(12)
#' n <- 300; z <- rnorm(n); v <- rnorm(n)
#' x <- z + v; u <- 0.5 * v + rnorm(n); y <- 1 + 0.5 * x + u
#' Zm <- cbind(1, z); Xm <- cbind(1, x)
#' beta <- solve(crossprod(Zm, Xm), crossprod(Zm, y))   # exactly identified IV
#' moments <- Zm * as.numeric(y - Xm %*% beta)          # z_t (y_t - x_t' beta)
#' D <- -crossprod(Zm, Xm) / n
#' fit <- aersn_gmm(moments, D, beta, target = 2, names = "slope")
#' aersn_test(fit, null = 0.5, draws = 2000)
#' @export
aersn_gmm <- function(moments, jacobian, estimate, weight = NULL,
                      target = NULL, profile = NULL, names = NULL) {
  M <- as_contribution_matrix(moments, "moments")
  n <- nrow(M)
  d <- ncol(M)
  D <- as.matrix(jacobian)
  storage.mode(D) <- "double"
  if (any(!is.finite(D))) aersn_abort("`jacobian` must be finite.")
  p <- ncol(D)
  if (p < 1L) aersn_abort("`jacobian` must have at least one column.")
  if (nrow(D) != d) {
    aersn_abort("`jacobian` has ", nrow(D), " rows but there are ", d,
                " moment conditions.")
  }
  if (d < p) {
    aersn_abort("Under-identified system: ", d, " moments for ", p,
                " parameters.")
  }
  if (qr(D)$rank < p) {
    aersn_abort("`jacobian` does not have full column rank (rank ",
                qr(D)$rank, " < p = ", p, "); the parameters are not ",
                "identified by these moments.")
  }
  beta <- as.numeric(estimate)
  if (length(beta) != p || any(!is.finite(beta))) {
    aersn_abort("`estimate` must be a finite vector of length p = ", p, ".")
  }
  if (!is.null(names(estimate))) names(beta) <- names(estimate)
  if (is.null(weight)) {
    W <- diag(d)
  } else {
    W <- as.matrix(weight)
    storage.mode(W) <- "double"
    if (!all(dim(W) == d)) aersn_abort("`weight` must be ", d, " x ", d, ".")
    if (any(!is.finite(W))) aersn_abort("`weight` must be finite.")
    if (max(abs(W - t(W))) > 1e-8 * max(1, max(abs(W)))) {
      aersn_abort("`weight` must be symmetric.")
    }
    if (min(eigen(W, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      aersn_abort("`weight` must be positive definite.")
    }
  }
  DW <- crossprod(D, W)                       # p x d
  Mhat <- solve(DW %*% D, DW)                 # p x d
  mbar <- colMeans(M)
  foc <- sqrt(n) * as.numeric(DW %*% mbar)
  scale_foc <- apply(M %*% t(DW), 2L, stats::sd)
  foc_standardized <- ifelse(scale_foc > 0, foc / scale_foc,
                             ifelse(foc == 0, 0, sign(foc) * Inf))
  foc_tolerance <- .Machine$double.eps^0.25
  if (any(abs(foc_standardized) > foc_tolerance)) {
    aersn_warn("The first-order condition is not close to zero relative ",
               "to the variation of the moment contributions (max ",
               "standardized residual = ",
               signif(max(abs(foc_standardized)), 3),
               "). Check optimizer convergence and that the estimate, ",
               "moments and weighting matrix correspond to the same GMM problem.")
  }
  psi_beta <- -M %*% t(Mhat)                  # rows: -M m_t
  tg <- resolve_target(target, beta, p, names)
  psi <- psi_beta %*% t(tg$J)
  aersn(tg$estimate, psi, profile = profile, names = tg$names,
        model = list(type = "gmm",
                     description = paste0("smooth GMM (", d, " moments, ", p,
                                          " parameters; psi_t = -h' M m_t)"),
                     foc = foc, foc_standardized = foc_standardized,
                     foc_tolerance = foc_tolerance,
                     p = p, d = d, target = tg$description),
        call = match.call())
}

## Interpret the `target` argument of the model interfaces.
resolve_target <- function(target, beta, p, names = NULL) {
  beta_names <- names(beta)
  if (is.null(beta_names)) beta_names <- paste0("beta", seq_len(p))
  if (is.null(target)) {
    J <- diag(p)
    est <- beta
    tnames <- beta_names
    desc <- "all parameters"
  } else if (is.list(target)) {
    if (!is.function(target$h)) {
      aersn_abort("A list `target` must contain a function `h`.")
    }
    est <- as.numeric(target$h(beta))
    J <- if (is.function(target$jacobian)) target$jacobian(beta) else
      target$jacobian
    if (is.null(J)) aersn_abort("A list `target` must contain `jacobian`.")
    J <- as.matrix(J)
    tnames <- if (!is.null(rownames(J))) rownames(J) else
      if (length(est) == 1L) "h(beta)" else paste0("h", seq_along(est))
    desc <- "smooth function h(beta)"
  } else if (is.matrix(target)) {
    J <- target
    storage.mode(J) <- "double"
    if (ncol(J) != p) aersn_abort("The target Jacobian must have ", p, " columns.")
    est <- as.numeric(J %*% beta)
    tnames <- if (!is.null(rownames(J))) rownames(J) else
      paste0("target", seq_len(nrow(J)))
    desc <- "linear transformation of beta"
  } else if (is.character(target)) {
    idx <- match(target, beta_names)
    if (anyNA(idx)) {
      aersn_abort("Unknown parameter name(s) in `target`: ",
                  paste(target[is.na(idx)], collapse = ", "))
    }
    idx <- check_indices(idx, p, "target")
    J <- diag(p)[idx, , drop = FALSE]
    est <- beta[idx]
    tnames <- beta_names[idx]
    desc <- paste0("selected parameters (", paste(tnames, collapse = ", "), ")")
  } else if (is.numeric(target)) {
    idx <- check_indices(target, p, "target")
    J <- diag(p)[idx, , drop = FALSE]
    est <- beta[idx]
    tnames <- beta_names[idx]
    desc <- paste0("selected parameters (", paste(tnames, collapse = ", "), ")")
  } else {
    aersn_abort("`target` must be NULL, indices, names, a matrix, or a list ",
                "with `h` and `jacobian`.")
  }
  if (any(!is.finite(J))) aersn_abort("The target Jacobian must be finite.")
  if (nrow(J) < 1L) aersn_abort("The target must have at least one coordinate.")
  if (ncol(J) != p) {
    aersn_abort("The target Jacobian has ", ncol(J), " columns but p = ", p)
  }
  if (length(est) != nrow(J)) {
    aersn_abort("The target has ", length(est), " values but its Jacobian ",
                "has ", nrow(J), " rows.")
  }
  if (qr(J)$rank < nrow(J)) {
    aersn_abort("The target transformation must have full row rank.")
  }
  if (!is.null(names)) {
    if (length(names) != nrow(J)) {
      aersn_abort("`names` must have length ", nrow(J), ".")
    }
    tnames <- names
  }
  list(estimate = stats::setNames(as.numeric(est), tnames), J = J,
       names = tnames, description = desc)
}

#' Linear regression by ordinary least squares
#'
#' Adjusted-range inference for regression coefficients.  OLS is the
#' exactly identified GMM estimator with moments
#' \eqn{m_t(b) = z_t(y_t - z_t^\top b)}; its influence contributions are
#' \eqn{\hat\psi_{n,t} = \hat{\dot h}_n(\hat Q^{-1} z_t \hat u_t)} with
#' \eqn{\hat Q = n^{-1}\sum_t z_t z_t^\top} and residuals \eqn{\hat u_t}
#' (manuscript Appendix B, OLS example).  The function forms the moments
#' and Jacobian from a fitted [stats::lm()] object and calls
#' [aersn_gmm()].
#'
#' @param object A formula (fitted with [stats::lm()] using `...`) or a
#'   fitted `lm` object.  Rows of the model frame must be in time order.
#' @param ... Passed to [stats::lm()] when `object` is a formula.
#' @param target Parameters of interest; see [aersn_gmm()].  The default is
#'   all coefficients including the intercept.
#' @param profile Optional variance-accumulation profile; see
#'   [aersn_gmm()].
#' @param names Optional names for the target parameters.
#' @return An [aersn()] object.
#' @details Weighted least squares, generalized linear models, multiple-response
#'   models and models with aliased (`NA`) coefficients are not supported.
#'   Observations removed by `na.action` (including `na.exclude`)
#'   are dropped before the path is formed; the function warns because the
#'   remaining rows are then treated as consecutive in time.  Validity
#'   requires the conditions of [aersn_gmm()] for the moments
#'   \eqn{z_t u_t}, in particular a functional central limit theorem for
#'   their partial sums and positive-definite \eqn{E(z_t z_t^\top)}.
#' @examples
#' set.seed(13)
#' n <- 250
#' x <- arima.sim(list(ar = 0.5), n)
#' u <- arima.sim(list(ar = 0.5), n)
#' y <- 1 + 0.8 * x + u
#' fit <- aersn_lm(y ~ x, target = "x")
#' aersn_test(fit, null = 0.8, draws = 2000)
#' fit2 <- aersn_lm(y ~ x)                           # both coefficients
#' confint(fit2, draws = 1000)
#' @export
aersn_lm <- function(object, ..., target = NULL, profile = NULL,
                     names = NULL) {
  if (inherits(object, "formula")) {
    fit <- stats::lm(object, ...)
  } else if (inherits(object, "lm")) {
    fit <- object
  } else {
    aersn_abort("`object` must be a formula or an lm fit.")
  }
  if (inherits(fit, "glm") || inherits(fit, "mlm")) {
    aersn_abort("`aersn_lm()` requires a single-response ordinary least-squares ",
                "fit; glm and mlm objects are not supported. Use influence ",
                "contributions or the score interface for other estimators.")
  }
  if (!is.null(fit$weights)) {
    aersn_abort("Weighted least squares is not supported in this version.")
  }
  b <- stats::coef(fit)
  if (anyNA(b)) {
    aersn_abort("The regression has aliased coefficients (",
                paste(names(b)[is.na(b)], collapse = ", "),
                "); remove collinear regressors.")
  }
  X <- stats::model.matrix(fit)
  ## The stored residuals follow the fitted rows; residuals(fit) pads them
  ## with NA under na.exclude, unlike model.matrix(fit).
  u <- fit$residuals
  if (length(u) != nrow(X)) {
    aersn_abort("The fitted residuals and design matrix have different row counts.")
  }
  if (!is.null(stats::na.action(fit))) {
    aersn_warn(length(stats::na.action(fit)), " observation(s) with missing ",
               "values were removed by lm(); the remaining rows are treated ",
               "as consecutive in time.")
  }
  n <- nrow(X)
  moments <- X * as.numeric(u)
  D <- -crossprod(X) / n
  out <- aersn_gmm(moments, D, b, weight = NULL, target = target,
                   profile = profile, names = names)
  out$model$type <- "lm"
  out$model$description <- paste0("linear regression by OLS (", n,
                                  " observations; psi_t = h' Q^{-1} z_t u_t)")
  out$model$formula <- stats::formula(fit)
  out$call <- match.call()
  out
}

#' Conditional likelihood scores
#'
#' Adjusted-range inference for a target \eqn{\theta = h(\beta)} of a
#' conditional maximum-likelihood (or M-) estimator from its fitted scores
#' and normalized negative Hessian.  The influence contributions are
#' \deqn{\hat\psi_{n,t} = \hat{\dot h}_n \hat J_n^{-1} s_t(\hat\beta_n),
#'   \qquad \hat J_n = n^{-1}J_{n,n}(\hat\beta_n),}
#' manuscript Appendix A.2 (equation (fitted-score-clock) and the text
#' following it).  Because the fitted scores sum to zero, the calendar-time
#' and profile-centered paths coincide; a profile changes only the
#' reference grid used for matched-grid quantiles.
#'
#' @details
#' With `profile = "opg"` the variance-accumulation profile is estimated
#' from cumulative score outer products by [aersn_opg_profile()].  This is
#' justified when the scores form a martingale-difference array and the
#' cumulative conditional score variance and cumulative expected negative
#' Hessian accumulate the same fraction \eqn{\tau(r)} of their full-sample
#' matrices (Supplement, Propositions "Multivariate fitted scores under a
#' common variance-accumulation profile" and "Feasible
#' variance-accumulation profile from score outer products").  Under
#' correct specification and information equality the two profiles agree.
#' For serially correlated scores the outer-product profile is not valid
#' and calendar-time centering should be used.  The proportionality measure
#' returned in `model$proportionality` describes departures from a common
#' scalar profile.
#'
#' @param scores `n x p` matrix of fitted conditional scores
#'   \eqn{s_t(\hat\beta_n)} in time order, one column per model parameter
#'   (including nuisance parameters).
#' @param neg_hessian `p x p` normalized negative Hessian
#'   \eqn{\hat J_n = -n^{-1}\sum_t \partial s_t(\hat\beta_n)/\partial\beta^\top}
#'   (or a consistent estimate of \eqn{J}); it must be symmetric and positive definite.
#' @param estimate The estimate \eqn{\hat\beta_n}.
#' @param target Parameters of interest; see [aersn_gmm()].
#' @param profile `"calendar"` (default), `"opg"` for the score
#'   outer-product profile, or a known profile as accepted by
#'   [aersn_profile()].
#' @param names Optional names for the target parameters.
#' @return An [aersn()] object.
#' @examples
#' ## Gaussian linear model with known variance: the score interface
#' ## reproduces the OLS influence contributions.
#' set.seed(14)
#' n <- 200; x <- rnorm(n); y <- 0.5 + x + rnorm(n)
#' X <- cbind(1, x); b <- solve(crossprod(X), crossprod(X, y))
#' s <- X * as.numeric(y - X %*% b)          # scores with sigma^2 = 1
#' J <- crossprod(X) / n
#' fit <- aersn_mle(s, J, b, target = 2, names = "slope")
#' confint(fit, draws = 4000)
#' fit_opg <- aersn_mle(s, J, b, target = 2, profile = "opg")
#' fit_opg$model$proportionality
#' @export
aersn_mle <- function(scores, neg_hessian, estimate, target = NULL,
                      profile = "calendar", names = NULL) {
  S <- as_contribution_matrix(scores, "scores")
  n <- nrow(S)
  p <- ncol(S)
  J <- as.matrix(neg_hessian)
  storage.mode(J) <- "double"
  if (!all(dim(J) == p)) aersn_abort("`neg_hessian` must be ", p, " x ", p)
  if (any(!is.finite(J))) aersn_abort("`neg_hessian` must be finite.")
  if (qr(J)$rank < p) aersn_abort("`neg_hessian` is singular.")
  if (max(abs(J - t(J))) > 1e-8 * max(abs(J))) {
    aersn_abort("`neg_hessian` must be symmetric.")
  }
  if (min(eigen(J, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    aersn_abort("`neg_hessian` must be positive definite.")
  }
  beta <- as.numeric(estimate)
  if (length(beta) != p || any(!is.finite(beta))) {
    aersn_abort("`estimate` must be a finite vector of length p = ", p, ".")
  }
  if (!is.null(names(estimate))) names(beta) <- names(estimate)
  score_sum <- colSums(S) / sqrt(n)
  if (max(abs(score_sum)) > 1e-6 * sqrt(n) * max(apply(S, 2L, stats::sd),
                                                 .Machine$double.eps)) {
    aersn_warn("The scaled full-sample score n^{-1/2} sum_t s_t is not ",
               "close to zero (max |.| = ", signif(max(abs(score_sum)), 3),
               "); the scores do not appear to be fitted scores at the ",
               "supplied estimate.")
  }
  prop <- NA_real_
  nodes <- NULL
  if (identical(profile, "calendar")) {
    nodes <- NULL
  } else if (identical(profile, "opg")) {
    nodes <- aersn_opg_profile(S)
    prop <- attr(nodes, "proportionality")
    nodes <- as.numeric(nodes)
  } else {
    nodes <- profile
  }
  psi_beta <- S %*% t(solve(J))               # rows: J^{-1} s_t
  tg <- resolve_target(target, beta, p, names)
  psi <- psi_beta %*% t(tg$J)
  aersn(tg$estimate, psi, profile = nodes, names = tg$names,
        model = list(type = "mle",
                     description = paste0("conditional likelihood scores (",
                                          p, " parameters; psi_t = h' J^{-1} s_t)",
                                          if (identical(profile, "opg"))
                                            "; outer-product profile" else ""),
                     score_sum = score_sum, proportionality = prop, p = p,
                     target = tg$description),
        call = match.call())
}
