#' Adjusted-range self-normalized inference from an estimate and its
#' influence contributions
#'
#' The central constructor of the package.  Given a parameter estimate and
#' the observation-level influence contributions of its asymptotically
#' linear representation, it builds the centered influence path, the
#' increment-hull self-normalizer and the diagnostics needed for tests
#' ([aersn_test()]), simultaneous confidence intervals ([confint.aersn()],
#' [aersn_contrast()]) and confidence regions ([aersn_region()]).  The
#' model interfaces [aersn_mean()], [aersn_lm()], [aersn_gmm()] and
#' [aersn_mle()] compute the contributions and call this function.
#'
#' @details
#' **Conventions.**  Parameters are ordered as in `estimate`; column `j`
#' of `psi` is the contribution to parameter `j`.  The representation
#' assumed is manuscript equation (asymptotic-linear-representation),
#' \deqn{\sqrt n(\hat\theta_n - \theta_0) = n^{-1/2}\sum_{t=1}^n \psi_t
#'   + o_p(1),}
#' with `n = nrow(psi)` observations in time order.  The centered path is
#' formed by [aersn_path()]: partial sums of the contributions scaled by
#' \eqn{n^{-1/2}} and centered by \eqn{\tau(k/n)} times the full-sample
#' sum, with \eqn{\tau(r) = r} under calendar-time centering.  The
#' statistic for a candidate value \eqn{v} is
#' \eqn{T_n(v) = \gamma_{K(\hat G_n)}\{\sqrt n(\hat\theta_n - v)\}} and the
#' level-\eqn{(1-\alpha)} confidence region is
#' \eqn{\hat\theta_n + (c_{q,1-\alpha}^{(n)}/\sqrt n) K(\hat G_n)}, where
#' \eqn{c_{q,1-\alpha}^{(n)}} is the matched-grid critical value from
#' [aersn_reference()].
#'
#' **Nuisance parameters and target transformations.**  When the
#' parameter of interest is a function \eqn{h(\beta)} of a larger estimated
#' vector, the contributions must already incorporate the first-order
#' effect of estimating the nuisance coordinates: for smooth GMM,
#' \eqn{\hat\psi_{n,t} = -\hat{\dot h}_n \hat M_n m_t(\hat\beta_n)}, and
#' for conditional likelihood
#' \eqn{\hat\psi_{n,t} = \hat{\dot h}_n \hat J_n^{-1} s_t(\hat\beta_n)}
#' (manuscript Appendix A).  Discarding nuisance coordinates from raw
#' moment contributions or scores does not give the correct path; use
#' [aersn_gmm()], [aersn_mle()] or [aersn_target()].
#'
#' **What the package does not check.**  Validity requires (i) a
#' functional central limit theorem for the partial sums of \eqn{\psi_t}
#' with a nonsingular scale matrix, (ii) the asymptotically linear
#' representation above, and (iii) Hausdorff convergence of the feasible
#' increment hull to that of the infeasible centered path (manuscript
#' Assumption 1).  Uniform convergence of the feasible path is sufficient
#' for (iii).  These conditions hold for the built-in interfaces under the
#' stated assumptions; for user-supplied contributions they are the user's
#' responsibility.  Under strong persistence the finite-sample null
#' rejection rate can exceed the nominal level (manuscript Section 5).
#'
#' @param estimate Numeric vector of length `q`: the estimate
#'   \eqn{\hat\theta_n} of the parameters of interest.
#' @param psi Numeric vector or `n x q` matrix of estimated influence
#'   contributions \eqn{\hat\psi_{n,t}}, observations in rows, in time
#'   order, columns ordered as `estimate`.
#' @param n Number of observations; must equal `nrow(psi)`.  It is the
#'   sample size of the asymptotically linear representation and enters the
#'   \eqn{n^{-1/2}} path scaling and the \eqn{\sqrt n} scaling of the
#'   estimation error.
#' @param profile `NULL` for calendar-time centering, or a
#'   variance-accumulation profile as accepted by [aersn_profile()].  See
#'   that help page for the additional conditions required.
#' @param names Optional parameter names.
#' @param model Optional list describing the estimator (used for printing);
#'   set by the model interfaces.
#' @param call The call to record.
#'
#' @return An object of class `"aersn"`: a list with elements `estimate`,
#'   `n`, `q`, `names`, `psi`, `path` (an [aersn_path()]), `hull` (an
#'   [aersn_hull()]), `centering`, `nodes`, `model` and `call`.
#' @seealso [aersn_test()], [confint.aersn()], [aersn_contrast()],
#'   [aersn_region()], [aersn_target()], [aersn_reference()].
#' @references Hong, Lin, Linton, Newey and Sun (2026), Sections 2-3 and
#'   Appendix A.
#' @examples
#' ## Synthetic example: mean of a bivariate AR(1) series
#' set.seed(10)
#' n <- 200; e <- matrix(rnorm(2 * n), n, 2); Y <- e
#' for (t in 2:n) Y[t, ] <- 0.4 * Y[t - 1, ] + e[t, ]
#' fit <- aersn(colMeans(Y), sweep(Y, 2, colMeans(Y)), names = c("m1", "m2"))
#' fit
#' aersn_test(fit, null = c(0, 0), draws = 500)          # small draws: example only
#' confint(fit, draws = 500)
#' @export
aersn <- function(estimate, psi, n = NROW(psi), profile = NULL, names = NULL,
                  model = NULL, call = NULL) {
  psi <- as_contribution_matrix(psi, "psi")
  if (!is.numeric(estimate) || any(!is.finite(estimate))) {
    aersn_abort("`estimate` must be a finite numeric vector.")
  }
  q <- ncol(psi)
  if (length(estimate) != q) {
    aersn_abort("`estimate` has length ", length(estimate), " but `psi` has ",
                q, " column(s); one column per parameter is required.")
  }
  check_scalar_number(n, "n", lower = 2, integer = TRUE)
  if (as.integer(n) != nrow(psi)) {
    aersn_abort("`n` = ", n, " but `psi` has ", nrow(psi), " rows. The ",
                "contributions must have one row per observation of the ",
                "asymptotically linear representation.")
  }
  names <- parameter_names(names, q, estimate = estimate, psi = psi)
  estimate <- stats::setNames(as.numeric(estimate), names)
  colnames(psi) <- names
  path <- aersn_path(psi, profile = profile, names = names)
  hull <- aersn_hull(path)
  zero <- which(hull$coordinate_ranges <= 0)
  if (length(zero)) {
    aersn_abort("The centered influence path has zero adjusted range for ",
                "parameter(s) ", paste(names[zero], collapse = ", "),
                ". This happens for a constant series or contributions ",
                "that are identically zero; no adjusted-range inference is ",
                "possible for these parameters.")
  }
  if (!hull$full_rank) {
    aersn_warn("The centered influence path is numerically rank deficient ",
               "(rank ", hull$rank, " < q = ", q, "). Joint statistics, ",
               "regions and simultaneous intervals cannot be evaluated; ",
               "marginal intervals remain available for coordinates with ",
               "positive adjusted range.")
  }
  if (is.null(model)) {
    model <- list(type = "influence",
                  description = "user-supplied influence contributions")
  }
  structure(list(estimate = estimate, n = as.integer(n), q = q,
                 names = names, psi = psi, path = path, hull = hull,
                 centering = path$centering, nodes = path$nodes,
                 model = model, call = if (is.null(call)) sys.call() else call),
            class = "aersn")
}

#' @export
print.aersn <- function(x, digits = 4L, ...) {
  cat("Affine-equivariant adjusted-range self-normalization (increment hull)\n")
  cat("  model:", x$model$description, "\n")
  cat("  n =", x$n, "observations; q =", x$q,
      if (x$q == 1L) "parameter\n" else "parameters\n")
  cat("  centering:", if (x$centering == "calendar")
    "calendar time, tau(r) = r" else
      "supplied variance-accumulation profile", "\n")
  cat("  estimate:\n")
  print(signif(x$estimate, digits))
  h <- x$hull
  cat("  hull diagnostics: numerical rank", h$rank, "of", x$q,
      if (h$full_rank) "" else "(RANK DEFICIENT)", "; condition",
      fmt_num(h$condition), "; min projected range",
      fmt_num(h$min_projected_range), "; spread",
      fmt_num(h$projected_range_spread), "\n")
  cat("Use aersn_test(), confint(), aersn_contrast(), aersn_region(), plot().\n")
  invisible(x)
}

#' @export
coef.aersn <- function(object, ...) object$estimate

#' @export
nobs.aersn <- function(object, ...) object$n

#' @rdname aersn_gauge
#' @export
aersn_gauge.aersn <- function(x, v, details = FALSE, ...) {
  V <- as_direction_matrix(v, x$q)
  Z <- sweep(-V, 2L, x$estimate, `+`) * sqrt(x$n)   # sqrt(n)(theta_hat - v)
  aersn_gauge(x$hull, Z, details = details)
}

## Signed scalar ratio z / range for q = 1 (used by one-sided tests).
signed_ratio <- function(object, v) {
  z <- sqrt(object$n) * (object$estimate - v)
  z / object$hull$coordinate_ranges[[1L]]
}

reference_for <- function(object, q, reference, draws, seed, cores = 1L,
                          verbose = FALSE) {
  resolve_reference(reference, q = q, n = object$n, nodes = object$nodes,
                    uniform = object$centering == "calendar",
                    draws = draws, seed = seed, cores = cores,
                    verbose = verbose)
}

#' Transform the target parameter
#'
#' Maps an [aersn()] object for a parameter vector \eqn{\beta} to one for a
#' target \eqn{\theta = h(\beta)} by transforming the estimate and the
#' influence contributions with the Jacobian of \eqn{h}:
#' \eqn{\hat\theta_n = h(\hat\beta_n)} and
#' \eqn{\hat\psi^\theta_{n,t} = \hat{\dot h}_n \hat\psi^\beta_{n,t}}.
#' For a linear transformation \eqn{A\beta}, the estimate is
#' \eqn{A\hat\beta_n} and the contributions are \eqn{A\hat\psi^\beta_{n,t}};
#' the increment hull transforms accordingly (\eqn{K(Ag) = AK(g)}).
#' Inference from the returned object uses the target dimension's reference
#' law. When dimension is reduced, its region generally differs from the
#' projection of the original region, which retains the original critical value.
#' Use [aersn_contrast()] with `type = "simultaneous"` for that projection.
#'
#' @param object An `"aersn"` object for the full parameter vector.
#' @param transform Either a numeric `m x q` matrix (linear target), or a
#'   function `h(beta)` returning the target vector.
#' @param jacobian For a function `transform`, its `m x q` Jacobian matrix
#'   at the estimate, or a function of `beta` returning it.
#' @param names Optional names for the target coordinates.
#' @return An `"aersn"` object for the target, with the same sample size and
#'   centering.
#' @details The transformation must have full row rank \eqn{m \le q};
#'   otherwise the transformed increment hull cannot be full dimensional.
#'   The delta method underlying the function-valued case requires `h` to
#'   be continuously differentiable at the true parameter with a
#'   full-row-rank derivative (manuscript Proposition "Feasible influence
#'   paths for smooth GMM").
#' @examples
#' set.seed(5)
#' fit <- aersn_mean(matrix(rnorm(300), 100, 3))
#' diff12 <- aersn_target(fit, matrix(c(1, -1, 0), 1, 3), names = "m1 - m2")
#' confint(diff12, draws = 500)
#' @export
aersn_target <- function(object, transform, jacobian = NULL, names = NULL) {
  check_aersn(object)
  beta <- object$estimate
  if (is.function(transform)) {
    est <- as.numeric(transform(beta))
    if (is.null(jacobian)) {
      aersn_abort("`jacobian` is required when `transform` is a function.")
    }
    J <- if (is.function(jacobian)) jacobian(beta) else jacobian
    J <- as.matrix(J)
  } else {
    J <- transform
    if (is.null(dim(J))) J <- matrix(as.numeric(J), nrow = 1L)
    J <- as.matrix(J)
    if (ncol(J) != object$q) {
      aersn_abort("The transformation has ", ncol(J), " columns but the ",
                  "object has q = ", object$q, " parameters.")
    }
    est <- as.numeric(J %*% beta)
  }
  storage.mode(J) <- "double"
  if (ncol(J) != object$q) {
    aersn_abort("The transformation has ", ncol(J), " columns but the ",
                "object has q = ", object$q, " parameters.")
  }
  m <- nrow(J)
  if (m < 1L) aersn_abort("The transformation must have at least one row.")
  if (length(est) != m) {
    aersn_abort("`transform` returned ", length(est), " values but its ",
                "Jacobian has ", m, " rows.")
  }
  if (any(!is.finite(J))) aersn_abort("The Jacobian must be finite.")
  if (qr(J)$rank < m) {
    aersn_abort("The transformation must have full row rank (rank ",
                qr(J)$rank, " < ", m, " rows); dependent contrasts cannot ",
                "define a full-dimensional target.")
  }
  if (is.null(names)) {
    names <- if (!is.null(rownames(J))) rownames(J) else
      if (m == 1L) "target" else paste0("target", seq_len(m))
  }
  psi_target <- object$psi %*% t(J)
  aersn(est, psi_target, n = object$n,
        profile = if (object$centering == "profile") object$nodes else NULL,
        names = names,
        model = list(type = "transformed",
                     description = paste0("transformed target (", m,
                                          " of ", object$q, " coordinates) of: ",
                                          object$model$description),
                     parent = object$model),
        call = match.call())
}

check_aersn <- function(object) {
  if (!inherits(object, "aersn")) {
    aersn_abort("`object` must be an aersn fit; see aersn().")
  }
  invisible(object)
}

require_full_rank <- function(object, what = "joint inference") {
  if (!object$hull$full_rank) {
    aersn_abort("The increment hull is not full dimensional (numerical rank ",
                object$hull$rank, " < q = ", object$q, "), so ", what,
                " is not evaluated.")
  }
  invisible(object)
}

#' Point-null test with the increment-hull gauge
#'
#' Tests \eqn{H_0: \theta_0 = v} (or \eqn{H_0: A\theta_0 = v} for a
#' full-row-rank contrast matrix \eqn{A}) using the statistic
#' \eqn{T_n(v) = \gamma_{K(\hat G_n)}\{\sqrt n(\hat\theta_n - v)\}} and the
#' reference law of [aersn_reference()].  The null is rejected at level
#' \eqn{\alpha = 1 -} `level` when \eqn{T_n(v)} exceeds the
#' \eqn{(1-\alpha)} quantile of the reference law.
#'
#' @details
#' A linear restriction \eqn{A\theta_0 = v} is tested as a point null for
#' the transformed target \eqn{A\theta} (see [aersn_target()]); nuisance
#' coordinates are handled through the influence contributions, not by
#' profiling.  Unrestricted composite hypotheses are not supported.  For
#' `q = 1` one-sided alternatives use the signed ratio
#' \eqn{\sqrt n(\hat\theta_n - v)} divided by the adjusted range, whose
#' limit is \eqn{M = Z/R} (Supplement, Section "Scalar specialization").
#'
#' The p-value is the reference tail probability at the observed statistic.
#' For a Monte Carlo two-sided reference it is the fraction of draws at or
#' above the statistic. One-sided scalar tails use the symmetry of the law,
#' consistently with the critical values. The result includes Monte Carlo
#' uncertainty; see [aersn_pvalue()]. For the closed-form scalar law no
#' Monte Carlo approximation is used. This is a limiting reference law,
#' not a finite-sample exact test for arbitrary dependent observations.
#'
#' @param object An [aersn()] object.
#' @param null Null value: a vector of length `q` (or of length `nrow(contrast)`).
#'   Defaults to zeros.
#' @param contrast Optional `m x q` matrix of full row rank defining the
#'   linear restriction `contrast %*% theta = null`.
#' @param level Confidence level \eqn{1 - \alpha} for the reported critical
#'   value.
#' @param alternative `"two.sided"`, or, for `q = 1` and `method = "hull"`,
#'   also `"greater"` or `"less"`.
#' @param method The inference method: a single identifier from
#'   [aersn_methods()] (`"hull"`, the default, `"ldl"`, `"shao"`, `"hac"`,
#'   `"fixedb"`, `"ewc"`), or a prebuilt [aersn_normalizer()] object.
#' @param ... Method-specific tuning arguments passed to
#'   [aersn_normalizer()], for example `kernel` and `bandwidth` for
#'   `method = "hac"`, `b` for `method = "fixedb"`, `nu` for
#'   `method = "ewc"`.  An unrecognized argument is an error rather than a
#'   silent default.
#' @param nu Number of cosine terms for `method = "ewc"`.  It has its own
#'   argument because `nu` is a prefix of `null` and would otherwise be
#'   captured by partial matching.
#' @param reference An [aersn_reference()] object matching the fit, the
#'   string `"continuous"` (closed-form law, `q = 1` only), the string
#'   `"registry"` (bundled tabulated quantiles, if available for this
#'   `q` and `n`), or `NULL` to simulate a matched-grid reference.
#' @param draws,seed,cores Settings passed to [aersn_reference()] when the
#'   reference is simulated. See its computational-cost discussion for
#'   memory requirements on long scalar series.
#' @return An object of class `c("aersn_test", "htest")` with the
#'   statistic, critical value, p-value, level, estimate, null value, the
#'   maximizing direction of the dual program, and the reference summary.
#' @examples
#' set.seed(6)
#' fit <- aersn_mean(rnorm(150, mean = 0.2))
#' aersn_test(fit, null = 0, reference = "continuous")
#' aersn_test(fit, null = 0, draws = 2000)
#' @export
aersn_test <- function(object, null = NULL, contrast = NULL, level = 0.95,
                       alternative = c("two.sided", "greater", "less"),
                       method = "hull", ..., nu = NULL,
                       reference = NULL, draws = NULL, seed = 1L,
                       cores = 1L) {
  check_aersn(object)
  alternative <- match.arg(alternative)
  tuning <- method_tuning(list(...), nu = nu, call = sys.call())
  check_level(level)
  data_name <- deparse1(substitute(object))
  if (!is.null(contrast)) {
    object <- aersn_target(object, contrast)
  }
  q <- object$q
  if (is.null(null)) null <- rep(0, q)
  null <- as.numeric(null)
  if (length(null) != q || any(!is.finite(null))) {
    aersn_abort("`null` must be a finite vector of length ", q, ".")
  }
  names(null) <- object$names
  nz <- resolve_method(object, method, tuning)
  require_method_full_rank(object, nz, "the joint test")
  if (alternative != "two.sided") {
    if (q != 1L) {
      aersn_abort("One-sided alternatives are available only for q = 1.")
    }
    if (nz$method != "hull") {
      aersn_abort("One-sided alternatives are available for method ",
                  "\"hull\" only. The reference law supplied for method '",
                  nz$method, "' is the law of its two-sided statistic; a ",
                  "signed law is not provided in this version.")
    }
  }
  z <- matrix(sqrt(object$n) * (object$estimate - null), nrow = 1L)
  stat <- normalizer_statistic(nz, z)
  direction <- if (nz$family == "geometric")
    aersn_gauge(nz$hull, z, details = TRUE)$direction else NULL
  ref <- reference_for_normalizer(nz, q, reference, object$n, object$nodes,
                                  object$centering == "calendar",
                                  draws, seed, cores)
  if (alternative == "two.sided") {
    cv <- aersn_critical_value(ref, level)
    pv <- aersn_pvalue(ref, stat)
    test_stat <- stat
  } else {
    if (level <= 0.5) aersn_abort("One-sided tests need level > 0.5.")
    signed <- signed_ratio(object, null)
    cv <- aersn_critical_value(ref, 2 * level - 1)
    if (alternative == "less") cv <- -cv
    pv <- aersn_pvalue(ref, signed, alternative = alternative)
    test_stat <- signed
  }
  title <- paste0(nz$label, " test, q = ", q)
  structure(list(statistic = c(T = test_stat),
                 parameter = c(q = q, n = object$n),
                 p.value = pv$p.value, critical.value = cv, level = level,
                 estimate = object$estimate, null.value = null,
                 alternative = alternative, method = title,
                 aersn_method = nz$method, tuning = nz$tuning,
                 statistic_scale = nz$scale, normalizer = nz,
                 data.name = data_name, direction = direction,
                 reference = reference_summary(ref), reference_type = ref$type,
                 p.mcse = pv$mcse, p.resolution = pv$resolution,
                 p.conf.int = pv$conf.int,
                 p.bracket = if (ref$type == "table") c(pv$lower, pv$upper) else NULL,
                 reject = if (alternative == "less") test_stat < cv else
                   test_stat > cv),
            class = c("aersn_test", "htest"))
}

#' @export
print.aersn_test <- function(x, digits = 4L, ...) {
  cat("\n\t", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat("T = ", fmt_num(x$statistic, digits), ", q = ", x$parameter[["q"]],
      ", n = ", x$parameter[["n"]], "\n", sep = "")
  cat("alternative hypothesis: ",
      switch(x$alternative,
             two.sided = "true parameter is not equal to the null value",
             greater = "true parameter is greater than the null value",
             less = "true parameter is less than the null value"), "\n",
      sep = "")
  cat("null value:", paste(names(x$null.value), "=",
                           fmt_num(x$null.value, digits), collapse = ", "),
      "\n")
  cat("estimate:  ", paste(names(x$estimate), "=",
                           fmt_num(x$estimate, digits), collapse = ", "),
      "\n")
  cat("critical value at level ", x$level, ": ", fmt_num(x$critical.value,
                                                          digits),
      if (x$reject) "  (reject)" else "  (do not reject)", "\n", sep = "")
  if (!is.null(x$tuning) && length(x$tuning)) {
    cat("tuning actually used: ", tuning_string(x$tuning), "\n", sep = "")
  }
  if (!is.null(x$p.bracket)) {
    cat(sprintf("p-value bracket from tabulated quantiles: %s < p <= %s\n",
                fmt_num(x$p.bracket[1L], 3), fmt_num(x$p.bracket[2L], 3)))
  } else if (x$reference_type == "monte_carlo") {
    if (x$p.value == 0) {
      cat("Estimated p-value = 0 (no reference draw exceeds the statistic)\n")
      if (!is.null(x$p.conf.int)) {
        cat(sprintf("95%% Monte Carlo interval for the reference tail probability: [%s, %s]\n",
                    fmt_num(x$p.conf.int[1L], 3), fmt_num(x$p.conf.int[2L], 3)))
      }
    } else {
      cat(sprintf("p-value = %s (Monte Carlo s.e. %s, resolution %s)\n",
                  fmt_num(x$p.value, 3), fmt_num(x$p.mcse, 2),
                  fmt_num(x$p.resolution, 2)))
    }
  } else {
    cat(sprintf("p-value = %s (closed-form law)\n", fmt_num(x$p.value, 4)))
  }
  cat("reference: ", x$reference, "\n", sep = "")
  invisible(x)
}

#' Confidence intervals from the increment hull
#'
#' Confidence intervals for coordinates of the parameter vector.  Three
#' constructions are distinguished:
#' * `"simultaneous"`: projections of the joint \eqn{q}-dimensional
#'   confidence region onto the coordinate axes, manuscript equation
#'   (projection-interval),
#'   \eqn{\hat\theta_j \pm (c_{q}/\sqrt n)\{\max_k \hat G_{n,k,j} -
#'   \min_k \hat G_{n,k,j}\}}, using the \eqn{q}-dimensional critical value.
#'   These intervals hold simultaneously for every linear projection of the
#'   region.
#' * `"joint"`: projections of the joint region of the selected
#'   coordinates only (dimension `m = length(parm)`), using \eqn{c_m}.
#'   Coverage is simultaneous over linear combinations of the selected
#'   coordinates.
#' * `"marginal"`: separately constructed scalar intervals, each applying
#'   the univariate adjusted-range statistic to the coordinate path with
#'   \eqn{c_1}.  They cover each coordinate individually and do not provide
#'   simultaneous coverage (manuscript Theorem 1(iii), last paragraph).
#'
#' @param object An [aersn()] object.
#' @param parm Coordinates to report: indices or names; default all.
#' @param level Confidence level.
#' @param type `"simultaneous"`, `"joint"` or `"marginal"`; see Details.
#' @param reference,draws,seed,cores Reference settings; see [aersn_test()].
#'   A supplied reference object must match the method, the tuning values and
#'   the dimension implied by `type`.
#' @param method The inference method: a single identifier from
#'   [aersn_methods()] (`"hull"`, the default, `"ldl"`, `"shao"`, `"hac"`,
#'   `"fixedb"`, `"ewc"`), or a prebuilt [aersn_normalizer()] object.
#' @param ... Method-specific tuning arguments passed to
#'   [aersn_normalizer()], for example `kernel` and `bandwidth` for
#'   `method = "hac"`, `b` for `method = "fixedb"`, `nu` for
#'   `method = "ewc"`.  An unrecognized argument is an error rather than a
#'   silent default.
#' @param nu Number of cosine terms for `method = "ewc"`.  It has its own
#'   argument because `nu` is a prefix of `null` and would otherwise be
#'   captured by partial matching.
#' @return A matrix with columns `lower` and `upper` and one row per
#'   coordinate, of class `"aersn_confint"`, with attributes `type`,
#'   `level`, `critical.value`, `dimension` (the dimension of the reference
#'   law used) and `reference`. For marginal inference, `tuning_by_contrast`
#'   records the tuning selected separately for each coordinate.
#' @examples
#' set.seed(7)
#' fit <- aersn_mean(matrix(rnorm(400), 200, 2))
#' confint(fit, draws = 1000)                      # simultaneous
#' confint(fit, type = "marginal", draws = 4000)   # separate scalar intervals
#' @export
confint.aersn <- function(object, parm = NULL, level = 0.95,
                          type = c("simultaneous", "joint", "marginal"),
                          method = "hull", ..., nu = NULL,
                          reference = NULL, draws = NULL, seed = 1L,
                          cores = 1L) {
  check_aersn(object)
  type <- match.arg(type)
  check_level(level)
  tuning <- method_tuning(list(...), nu = nu, call = sys.call())
  idx <- resolve_parm(parm, object$names)
  A <- selection_matrix(idx, object$q)
  parts <- method_intervals(object, A, object$names[idx], type, method,
                            tuning, level, reference, draws, seed, cores)
  est <- object$estimate[idx]
  out <- cbind(lower = est - parts$half, upper = est + parts$half)
  rownames(out) <- object$names[idx]
  structure(out, type = type, level = level, critical.value = parts$cv,
            dimension = parts$dim, reference = reference_summary(parts$ref),
            aersn_method = parts$nz$method, tuning = parts$nz$tuning,
            tuning_by_contrast = parts$tuning_by_contrast,
            class = c("aersn_confint", "matrix", "array"))
}

resolve_parm <- function(parm, names) {
  if (is.null(parm)) return(seq_along(names))
  if (is.character(parm)) {
    idx <- match(parm, names)
    if (anyNA(idx)) {
      aersn_abort("Unknown parameter name(s): ",
                  paste(parm[is.na(idx)], collapse = ", "))
    }
    return(check_indices(idx, length(names), "parm"))
  }
  check_indices(parm, length(names), "parm")
}

#' @export
print.aersn_confint <- function(x, digits = 4L, ...) {
  type <- attr(x, "type")
  cat(sprintf("%s %s confidence intervals, level %s\n",
              switch(type, simultaneous = "Simultaneous (joint-region projection)",
                     joint = "Joint-region projection (selected coordinates)",
                     marginal = "Marginal (separately constructed)"),
              result_method_label(x), attr(x, "level")))
  plain <- matrix(as.numeric(x), nrow = nrow(x), dimnames = dimnames(x))
  print(signif(plain, digits))
  cat(sprintf("critical value %s (reference dimension %d); %s\n",
              fmt_num(attr(x, "critical.value"), digits),
              attr(x, "dimension"), attr(x, "reference")))
  if (type == "marginal") {
    cat("Note: marginal intervals do not provide simultaneous coverage.\n")
  }
  invisible(x)
}

#' Simultaneous intervals for linear contrasts
#'
#' Confidence intervals for the contrasts \eqn{a_k^\top\theta}, rows of a
#' contrast matrix, obtained from the support function of the confidence
#' region: \eqn{a_k^\top\hat\theta_n \pm (c/\sqrt n)\,h_{K(\hat G_n)}(a_k)},
#' where \eqn{h_K(a)} is the adjusted range of the projected path
#' \eqn{a^\top\hat G_n} (manuscript equation (projection-interval)).
#'
#' @inheritParams confint.aersn
#' @param contrasts A vector of length `q` or an `m x q` matrix with one
#'   contrast per row.  Row names are used as labels.
#' @param type `"simultaneous"` uses the \eqn{q}-dimensional critical value
#'   (projections of the full joint region; simultaneous over all linear
#'   contrasts); `"joint"` uses the \eqn{m}-dimensional critical value
#'   (projections of the joint region of the `m` specified contrasts, which
#'   must have full row rank); `"marginal"` uses the scalar critical value
#'   for each contrast separately.
#' @return A data frame with columns `estimate`, `lower`, `upper`, one row
#'   per contrast, with attributes `type`, `level`, `critical.value`,
#'   `dimension` and `reference`. For marginal inference, `tuning_by_contrast`
#'   records the tuning selected separately for each contrast.
#' @examples
#' set.seed(8)
#' fit <- aersn_mean(matrix(rnorm(300), 100, 3), names = c("a", "b", "c"))
#' A <- rbind("a - b" = c(1, -1, 0), "mean" = c(1, 1, 1) / 3)
#' aersn_contrast(fit, A, draws = 500)
#' @export
aersn_contrast <- function(object, contrasts, level = 0.95,
                           type = c("simultaneous", "joint", "marginal"),
                           method = "hull", ..., nu = NULL,
                           reference = NULL, draws = NULL, seed = 1L,
                           cores = 1L) {
  check_aersn(object)
  type <- match.arg(type)
  check_level(level)
  tuning <- method_tuning(list(...), nu = nu, call = sys.call())
  A <- as_direction_matrix(contrasts, object$q)
  m <- nrow(A)
  labels <- if (!is.null(rownames(A))) rownames(A) else
    paste0("contrast", seq_len(m))
  parts <- method_intervals(object, A, labels, type, method, tuning, level,
                            reference, draws, seed, cores)
  est <- parts$estimate
  out <- data.frame(estimate = est, lower = est - parts$half,
                    upper = est + parts$half, row.names = labels)
  structure(out, type = type, level = level, critical.value = parts$cv,
            dimension = parts$dim, reference = reference_summary(parts$ref),
            aersn_method = parts$nz$method, tuning = parts$nz$tuning,
            tuning_by_contrast = parts$tuning_by_contrast,
            class = c("aersn_contrast", "data.frame"))
}

## Subsetting a result table keeps its class only while the table still has
## the columns its print method needs. Dropping columns returns a plain data
## frame, which prints correctly, instead of a classed object whose
## attributes no longer describe it.
#' @export
`[.aersn_contrast` <- function(x, ...) {
  drop_result_class(NextMethod(), x, "aersn_contrast",
                    c("type", "level", "critical.value", "dimension",
                      "reference", "aersn_method", "tuning", "tuning_by_contrast"))
}

drop_result_class <- function(out, x, cls, attrs) {
  if (!inherits(out, cls)) return(out)
  if (identical(names(out), names(x))) return(out)
  class(out) <- setdiff(class(out), cls)
  for (a in attrs) attr(out, a) <- NULL
  out
}

#' @export
print.aersn_contrast <- function(x, digits = 4L, ...) {
  if (!all(c("estimate", "lower", "upper") %in% names(x))) {
    return(print(as.data.frame(unclass(x)), digits = digits))
  }
  type <- attr(x, "type")
  cat(sprintf("%s %s intervals for linear contrasts, level %s\n",
              switch(type, simultaneous = "Simultaneous (joint-region projection)",
                     joint = "Joint-region projection (specified contrasts)",
                     marginal = "Marginal (separately constructed)"),
              result_method_label(x), attr(x, "level")))
  print(as.data.frame(lapply(unclass(x), signif, digits = digits),
                      row.names = row.names(x)))
  cat(sprintf("critical value %s (reference dimension %d); %s\n",
              fmt_num(attr(x, "critical.value"), digits),
              attr(x, "dimension"), attr(x, "reference")))
  if (type == "marginal") {
    cat("Note: marginal intervals do not provide simultaneous coverage.\n")
  }
  invisible(x)
}

#' Joint confidence region
#'
#' The level-\eqn{(1-\alpha)} confidence region
#' \eqn{\mathcal C_n = \hat\theta_n + (c^{(n)}_{q,1-\alpha}/\sqrt n)
#' K(\hat G_n)} (manuscript equation (confidence-region)).  The region is a
#' translated and scaled copy of the increment hull; it is convex,
#' centrally symmetric about the estimate, generally not an ellipsoid, and
#' affine equivariant.  It is represented by its center, scale, and hull;
#' membership is evaluated through the gauge and projections through the
#' support function.  Explicit vertices are available for `q = 2`
#' ([aersn_vertices()]), and two-dimensional projections or slices for
#' plotting in any dimension ([plot.aersn_region()]).
#'
#' @inheritParams confint.aersn
#' @return An object of class `"aersn_region"` with elements `center`,
#'   `scale` (\eqn{c/\sqrt n}), `critical.value`, `level`, `hull`, `q`,
#'   `n`, `names` and `reference`.
#' @seealso [aersn_contains()], [aersn_support()], [aersn_gauge()],
#'   [aersn_vertices()], [plot.aersn_region()].
#' @examples
#' set.seed(9)
#' fit <- aersn_mean(matrix(rnorm(400), 200, 2))
#' reg <- aersn_region(fit, draws = 1000)
#' reg
#' aersn_contains(reg, rbind(c(0, 0), c(0.5, 0.5)))
#' plot(reg)
#' @export
aersn_region <- function(object, level = 0.95, method = "hull", ..., nu = NULL,
                         reference = NULL, draws = NULL, seed = 1L,
                         cores = 1L) {
  check_aersn(object)
  check_level(level)
  nz <- resolve_method(object, method,
                       method_tuning(list(...), nu = nu, call = sys.call()))
  require_method_full_rank(object, nz, "the confidence region")
  ref <- reference_for_normalizer(nz, object$q, reference, object$n,
                                  object$nodes,
                                  object$centering == "calendar",
                                  draws, seed, cores)
  cv <- aersn_critical_value(ref, level)
  structure(list(center = object$estimate, scale = cv / sqrt(object$n),
                 critical.value = cv, level = level,
                 hull = if (nz$family == "geometric") nz$hull else NULL,
                 normalizer = nz, method = nz$method, family = nz$family,
                 q = object$q, n = object$n, names = object$names,
                 centering = object$centering,
                 reference = reference_summary(ref), reference_object = ref),
            class = "aersn_region")
}

#' @export
print.aersn_region <- function(x, digits = 4L, ...) {
  cat(x$normalizer$label, " confidence region, level ", x$level, "\n", sep = "")
  cat("  q =", x$q, " n =", x$n, " critical value", fmt_num(x$critical.value,
                                                             digits),
      " scale c/sqrt(n) =", fmt_num(x$scale, digits), "\n")
  cat("  center (estimate):\n")
  print(signif(x$center, digits))
  cat("  coordinate projections (simultaneous):\n")
  ci <- unclass(confint(x))
  attributes(ci) <- attributes(ci)[c("dim", "dimnames")]
  print(signif(ci, digits))
  if (x$q == 2L) {
    if (x$family == "geometric") {
      V <- aersn_vertices(x)
      cat("  area:", fmt_num(polygon_area(V), digits), "with", nrow(V),
          "vertices (convex polygon)\n")
    } else {
      cat("  area:", fmt_num(region_area_2d(x), digits), "(ellipse)\n")
    }
  }
  if (length(x$normalizer$tuning)) {
    cat("  tuning actually used:", tuning_string(x$normalizer$tuning), "\n")
  }
  cat("  reference:", x$reference, "\n")
  invisible(x)
}

## Area of a two-dimensional region.
region_area_2d <- function(region) {
  if (region$q != 2L) aersn_abort("An area is reported only for q = 2.")
  if (region$family == "geometric") {
    return(polygon_area(aersn_vertices(region)))
  }
  V <- region$normalizer$matrix
  pi * (region$critical.value / region$n) * sqrt(det(V))
}

## Gauge of the region itself: the factor by which the region must be scaled
## about its centre to contain v.  It equals one on the boundary for both
## families, and is homogeneous of degree one in v - centre.
region_gauge <- function(region, v) {
  V <- as_direction_matrix(v, region$q)
  z <- sweep(V, 2L, region$center, `-`) * sqrt(region$n)
  stat <- normalizer_statistic(region$normalizer, z)
  if (region$normalizer$scale == "gauge") {
    stat / region$critical.value
  } else {
    sqrt(pmax(stat, 0) / region$critical.value)
  }
}

#' @rdname aersn_gauge
#' @export
aersn_gauge.aersn_region <- function(x, v, details = FALSE, ...) {
  if (x$family == "geometric") {
    V <- as_direction_matrix(v, x$q)
    X <- sweep(V, 2L, x$center, `-`) / x$scale
    return(aersn_gauge(x$hull, X, details = details))
  }
  if (details) {
    aersn_abort("`details = TRUE` reports the maximizing direction of the ",
                "increment-hull linear program and applies to method ",
                "\"hull\" only.")
  }
  region_gauge(x, v)
}

#' @rdname aersn_support
#' @export
aersn_support.aersn_region <- function(x, u, ...) {
  U <- as_direction_matrix(u, x$q)
  as.numeric(U %*% x$center) +
    normalizer_radius(x$normalizer, U, x$critical.value)
}

#' Membership test for a confidence region
#'
#' Tests whether candidate parameter vectors belong to a confidence region,
#' equivalently whether the point-null test does not reject them.
#'
#' @param region An [aersn_region()] object.
#' @param v A vector of length `q` or a matrix with one candidate per row.
#' @param tol Numerical tolerance on the gauge.
#' @return A logical vector with attribute `"gauge"` holding the gauge of
#'   each candidate relative to the region (at most one inside).
#' @export
aersn_contains <- function(region, v, tol = 1e-9) {
  check_scalar_number(tol, "tol", lower = 0)
  if (!inherits(region, "aersn_region")) {
    aersn_abort("`region` must be an aersn_region object.")
  }
  g <- aersn_gauge(region, v)
  structure(g <= 1 + tol, gauge = g)
}

#' @export
confint.aersn_region <- function(object, parm = NULL, level = NULL, ...) {
  if (!is.null(level) && !isTRUE(all.equal(level, object$level))) {
    aersn_abort("The region was built at level ", object$level,
                "; rebuild it with aersn_region() for another level.")
  }
  idx <- resolve_parm(parm, object$names)
  half <- normalizer_radius(object$normalizer,
                            selection_matrix(idx, object$q),
                            object$critical.value)
  est <- object$center[idx]
  out <- cbind(lower = est - half, upper = est + half)
  rownames(out) <- object$names[idx]
  structure(out, type = "simultaneous", level = object$level,
            critical.value = object$critical.value, dimension = object$q,
            reference = object$reference,
            aersn_method = object$method,
            tuning = object$normalizer$tuning,
            class = c("aersn_confint", "matrix", "array"))
}

#' Vertices of a two-dimensional confidence region
#'
#' For `q = 2`, the increment hull is a polygon whose vertices are among the
#' pairwise differences of the vertices of the convex hull of the path
#' values (\eqn{K(g) = P(g) - P(g)}, Supplement Proposition
#' "Characterization of the bridge-increment hull").  The region's vertices
#' are the center plus the scaled hull vertices, in counter-clockwise order.
#'
#' @param region An [aersn_region()] object with `q = 2`, or an
#'   [aersn_hull()] object with `q = 2` (then the unscaled hull vertices are
#'   returned).
#' @return A numeric matrix with two columns.
#' @export
aersn_vertices <- function(region) {
  if (inherits(region, "aersn_hull")) {
    if (region$q != 2L) aersn_abort("Vertices are available only for q = 2.")
    return(hull_vertices_2d(region$G))
  }
  if (!inherits(region, "aersn_region")) {
    aersn_abort("`region` must be an aersn_region or aersn_hull object.")
  }
  if (region$q != 2L) {
    aersn_abort("Vertices are available only for q = 2; use ",
                "aersn_projection() or aersn_slice() for two-dimensional ",
                "views of higher-dimensional regions.")
  }
  if (!identical(region$family, "geometric")) {
    aersn_abort("Region vertices exist only for the increment hull, whose ",
                "region is a convex polygon. The region for method '",
                region$method, "' is an ellipsoid; use aersn_projection() ",
                "or aersn_slice() to obtain boundary points.")
  }
  V <- hull_vertices_2d(region$hull$G)
  sweep(V * region$scale, 2L, region$center, `+`)
}

polygon_area <- function(V) {
  m <- nrow(V)
  nxt <- c(seq_len(m)[-1L], 1L)
  abs(sum(V[, 1L] * V[nxt, 2L] - V[nxt, 1L] * V[, 2L])) / 2
}

#' Two-dimensional projections and slices of a confidence region
#'
#' `aersn_projection()` returns the polygon of the orthogonal projection of
#' the region onto two coordinates.  Because projection is linear, the
#' projection of \eqn{\hat\theta_n + s K(\hat G_n)} is the same construction
#' applied to the two-dimensional sub-path, so it is exact and has
#' simultaneous coverage.  `aersn_slice()` returns the boundary of the
#' intersection of the region with the plane through the center in which
#' all other coordinates are held at their estimates; each boundary point is
#' \eqn{\hat\theta_n + s\,d/\gamma_K(d)} for a direction \eqn{d} in the
#' plane, evaluated by the gauge linear program.  A slice is a
#' cross-section, not a confidence set for the two coordinates.
#'
#' @param region An [aersn_region()] object.
#' @param coords Integer vector of length two: the coordinates to display.
#' @param n_angles Number of directions used for the slice boundary.
#' @return A numeric matrix with two columns (polygon vertices in order).
#' @export
aersn_projection <- function(region, coords = c(1L, 2L)) {
  if (!inherits(region, "aersn_region")) {
    aersn_abort("`region` must be an aersn_region object.")
  }
  coords <- check_indices(coords, region$q, "coords")
  if (length(coords) != 2L || any(coords < 1L | coords > region$q) ||
      coords[1L] == coords[2L]) {
    aersn_abort("`coords` must be two distinct coordinates in 1..", region$q)
  }
  if (region$family == "geometric") {
    G2 <- region$hull$G[, coords, drop = FALSE]
    V <- hull_vertices_2d(G2)
    return(sweep(V * region$scale, 2L, region$center[coords], `+`))
  }
  ## The projection of an ellipsoid onto two coordinates is the ellipse of
  ## the corresponding submatrix of the normalizer.
  V22 <- region$normalizer$matrix[coords, coords, drop = FALSE]
  L <- t(chol((V22 + t(V22)) / 2))
  phi <- seq(0, 2 * pi, length.out = 181L)[-181L]
  B <- cbind(cos(phi), sin(phi)) %*% t(L) *
    sqrt(region$critical.value / region$n)
  out <- sweep(B, 2L, region$center[coords], `+`)
  colnames(out) <- region$names[coords]
  out
}

#' @rdname aersn_projection
#' @export
aersn_slice <- function(region, coords = c(1L, 2L), n_angles = 180L) {
  if (!inherits(region, "aersn_region")) {
    aersn_abort("`region` must be an aersn_region object.")
  }
  coords <- check_indices(coords, region$q, "coords")
  if (length(coords) != 2L || any(coords < 1L | coords > region$q) ||
      coords[1L] == coords[2L]) {
    aersn_abort("`coords` must be two distinct coordinates in 1..", region$q)
  }
  check_scalar_number(n_angles, "n_angles", lower = 3, integer = TRUE)
  phi <- seq(0, 2 * pi, length.out = n_angles + 1L)[-(n_angles + 1L)]
  D <- matrix(0, n_angles, region$q)
  D[, coords[1L]] <- cos(phi)
  D[, coords[2L]] <- sin(phi)
  ## The boundary in direction d is centre + d / gamma(centre + d), where
  ## gamma is the region's own gauge; this holds for both families.
  g <- region_gauge(region, sweep(D, 2L, region$center, `+`))
  B <- D[, coords, drop = FALSE] / g
  out <- sweep(B, 2L, region$center[coords], `+`)
  colnames(out) <- region$names[coords]
  out
}

#' @export
summary.aersn <- function(object, level = 0.95, method = "hull", ...,
                          reference = NULL, draws = NULL, seed = 1L,
                          cores = 1L) {
  check_aersn(object)
  ## The rank requirement differs by method, so let confint() decide and
  ## report its reason rather than gating everything on the hull's rank.
  ci <- tryCatch(
    confint(object, level = level, method = method, ...,
            reference = reference, draws = draws, seed = seed, cores = cores),
    error = function(e) structure(character(), reason = conditionMessage(e)))
  structure(list(fit = object, level = level, confint = ci),
            class = "summary.aersn")
}

#' @export
print.summary.aersn <- function(x, digits = 4L, ...) {
  print(x$fit, digits = digits)
  if (inherits(x$confint, "aersn_confint")) {
    cat("\n")
    print(x$confint, digits = digits)
  } else if (!is.null(attr(x$confint, "reason"))) {
    cat("\nNo confidence intervals: ", attr(x$confint, "reason"), "\n",
        sep = "")
  }
  invisible(x)
}


result_method_label <- function(x) {
  method <- attr(x, "aersn_method")
  if (is.null(method)) method <- "hull"
  method_spec(method)$label
}
