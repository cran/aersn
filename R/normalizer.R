## Normalizer layer shared by every inference method.
##
## A normalizer turns the scaled estimation error z = sqrt(n)(theta_hat - v)
## into a test statistic, and determines the shape of the confidence region.
## Two families are supported.
##
##   family = "geometric": the increment hull of the centered influence path.
##     The statistic is the gauge gamma_K(z), homogeneous of degree one in z.
##     The region theta_hat + (c/sqrt n) K is generally not an ellipsoid.
##
##   family = "quadratic": a q x q positive-definite matrix V.
##     The statistic is the Wald form z' V^{-1} z, homogeneous of degree two.
##     The region is the ellipsoid {v : n (theta_hat - v)' V^{-1} (theta_hat - v) <= c}.
##
## Because the two families have different homogeneity, their statistics and
## critical values are on different scales and must never be compared as
## numbers; only the resulting tests, p-values and interval widths are
## comparable.  Every normalizer carries the reference family it requires and
## the tuning values actually used, so that a statistic can never be combined
## with a reference law computed for a different method or tuning value.

## Registry of supported methods.  `needs_calendar` marks methods whose
## discretization is defined on the calendar-time grid.
aersn_method_table <- function() {
  list(
    hull = list(
      label = "Adjusted-range increment hull",
      family = "geometric", scale = "gauge",
      reference = "hull_gauge", needs_calendar = FALSE,
      self_normalized = TRUE),
    ldl = list(
      label = "LDL componentwise adjusted range",
      family = "quadratic", scale = "wald",
      reference = "componentwise_mq2", needs_calendar = FALSE,
      self_normalized = TRUE),
    shao = list(
      label = "Quadratic self-normalization (Shao)",
      family = "quadratic", scale = "wald",
      reference = "shao_sq2", needs_calendar = FALSE,
      self_normalized = TRUE),
    hac = list(
      label = "HAC long-run covariance",
      family = "quadratic", scale = "wald",
      reference = "chisq", needs_calendar = TRUE,
      self_normalized = FALSE),
    fixedb = list(
      label = "Bartlett fixed-b",
      family = "quadratic", scale = "wald",
      reference = "fixedb_bartlett", needs_calendar = TRUE,
      self_normalized = FALSE),
    ewc = list(
      label = "Equal-weighted cosine (EWC)",
      family = "quadratic", scale = "wald",
      reference = "ewc_F", needs_calendar = TRUE,
      self_normalized = FALSE)
  )
}

#' Method identifiers
#'
#' The character identifiers accepted by the `method` argument of
#' [aersn_test()], [confint.aersn()], [aersn_contrast()], [aersn_region()]
#' and [aersn_compare()], in the order used by comparison tables.
#'
#' @details
#' * `"hull"`: the affine-equivariant adjusted-range increment hull, the
#'   method of the package's main reference and the default everywhere.
#' * `"ldl"`: componentwise adjusted ranges after lag-zero partial
#'   prewhitening by an LDL factorization ([aersn_ldl_normalizer()]).
#' * `"shao"`: quadratic self-normalization ([aersn_shao_normalizer()]).
#' * `"hac"`: kernel long-run covariance estimation ([aersn_hac_lrv()]).
#' * `"fixedb"`: Bartlett fixed-b ([aersn_fixed_b_normalizer()]).
#' * `"ewc"`: equal-weighted cosine ([aersn_ewc_lrv()]).
#'
#' `"hull"`, `"ldl"` and `"shao"` are self-normalized: the normalizer is a
#' random object with a nondegenerate limit and the reference law is
#' nonstandard.  `"hac"` uses a consistent covariance estimate with a
#' chi-squared reference.  `"fixedb"` and `"ewc"` use inconsistent covariance
#' estimates with fixed-smoothing reference laws.
#'
#' @return A character vector of method identifiers.
#' @examples
#' aersn_methods()
#' @export
aersn_methods <- function() names(aersn_method_table())

method_spec <- function(method) {
  tab <- aersn_method_table()
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    aersn_abort("`method` must be a single method identifier, one of ",
                paste(names(tab), collapse = ", "), ".")
  }
  if (!method %in% names(tab)) {
    aersn_abort("Unknown method '", method, "'. Supported methods are ",
                paste(names(tab), collapse = ", "), ".")
  }
  tab[[method]]
}

## Build a normalizer object.  Internal: the exported constructors below call
## this after computing their method-specific parts.
new_normalizer <- function(method, q, n, names, matrix = NULL, hull = NULL,
                           tuning = list(), diagnostics = list(),
                           reference_args = list(), centering = "calendar",
                           nodes = NULL, notes = character(), source) {
  spec <- method_spec(method)
  check_method_grid(method, source)
  structure(list(source = list(psi = source$psi, nodes = source$nodes,
                               centering = source$centering),
                 method = method, label = spec$label, family = spec$family,
                 scale = spec$scale, reference_family = spec$reference,
                 self_normalized = spec$self_normalized,
                 q = q, n = n, names = names,
                 matrix = matrix, hull = hull,
                 tuning = tuning, diagnostics = diagnostics,
                 reference_args = reference_args,
                 centering = centering, nodes = nodes, notes = notes),
            class = c(paste0("aersn_normalizer_", method), "aersn_normalizer"))
}

#' Construct the normalizer used by an inference method
#'
#' Builds the object that a method uses to normalize the estimation error:
#' the increment hull for `"hull"`, or a `q` by `q` positive-definite matrix
#' for the other five methods.  The inference functions call this internally;
#' it is exported so that a normalizer can be inspected, reused across
#' several tests, or built once and passed to [aersn_test()] and friends.
#' A prebuilt normalizer is tied to its influence contributions and profile.
#' For a different data set or a lower-dimensional target, pass a method name
#' and tuning arguments to construct the normalizer for that target.
#'
#' @details
#' Method-specific tuning arguments are passed through `...` and validated by
#' the individual constructors:
#' * `"hull"`: none.
#' * `"ldl"`: `order`.
#' * `"shao"`: `integration`.
#' * `"hac"`: `kernel`, `bandwidth`, `lag`, `prewhite`, `adjust`, `center`.
#' * `"fixedb"`: `b`.
#' * `"ewc"`: `nu`, `center`.
#'
#' Passing an unknown tuning argument is an error rather than a silent
#' default, so that a misspelled kernel or bandwidth rule cannot be ignored.
#'
#' @param object An [aersn()] fit.
#' @param method A method identifier; see [aersn_methods()].
#' @param ... Method-specific tuning arguments.
#'
#' @return An object of class `"aersn_normalizer"`.  Common fields are
#'   `method`, `label`, `family` (`"geometric"` or `"quadratic"`), `scale`
#'   (`"gauge"` or `"wald"`), `reference_family`, `q`, `n`, `names`,
#'   `matrix` (the long-run covariance or quadratic normalizer, `NULL` for
#'   the hull), `hull` (for the hull method), `tuning` (the realized tuning
#'   values) and `diagnostics`.
#' @seealso [aersn_ldl_normalizer()], [aersn_shao_normalizer()],
#'   [aersn_hac_lrv()], [aersn_fixed_b_normalizer()], [aersn_ewc_lrv()],
#'   [aersn_compare()].
#' @examples
#' set.seed(1)
#' fit <- aersn_mean(matrix(rnorm(300), 150, 2))
#' aersn_normalizer(fit, "shao")
#' aersn_normalizer(fit, "hac", kernel = "Parzen", bandwidth = "andrews")
#' @export
aersn_normalizer <- function(object, method = "hull", ...) {
  check_aersn(object)
  if (inherits(method, "aersn_normalizer")) {
    if (length(list(...))) {
      aersn_abort("Tuning arguments cannot be combined with a prebuilt normalizer.")
    }
    check_normalizer_fit(method, object)
    return(method)
  }
  spec <- method_spec(method)
  if (spec$needs_calendar && object$centering != "calendar") {
    aersn_abort("Method '", method, "' is defined on the calendar-time grid, ",
                "but this fit uses a supplied variance-accumulation profile. ",
                "Use method = \"hull\", \"ldl\", or \"shao\" with ",
                "integration = \"profile\" for profile-centered fits.")
  }
  nz <- switch(method,
               hull = aersn_hull_normalizer(object, ...),
               ldl = aersn_ldl_normalizer(object, ...),
               shao = aersn_shao_normalizer(object, ...),
               hac = aersn_hac_lrv(object, ...),
               fixedb = aersn_fixed_b_normalizer(object, ...),
               ewc = aersn_ewc_lrv(object, ...))
  ## Record the arguments as supplied, so that the same method with the same
  ## tuning can be rebuilt on a lower-dimensional target.
  nz$call_args <- list(...)
  nz
}

## Rebuild a normalizer with the same method and the same supplied tuning on
## another fit, typically a linear target of the original one.  Data-dependent
## tuning, such as an automatic HAC bandwidth, is re-selected from the
## target's own contributions; the realized values are recorded again.
rebuild_normalizer <- function(nz, object) {
  args <- if (is.null(nz$call_args)) list() else nz$call_args
  do.call(aersn_normalizer, c(list(object, nz$method), args))
}

## Guard against combining a normalizer with a different fit.
check_normalizer_fit <- function(nz, object) {
  if (nz$q != object$q || nz$n != object$n) {
    aersn_abort("The supplied normalizer was built for q = ", nz$q, " and n = ",
                nz$n, ", but the fit has q = ", object$q, " and n = ",
                object$n, ".")
  }
  if (!identical(nz$names, object$names)) {
    aersn_abort("The supplied normalizer was built for parameters ",
                paste(nz$names, collapse = ", "), " but the fit has ",
                paste(object$names, collapse = ", "), ".")
  }
  if (is.null(nz$source) || !identical(nz$source$psi, object$psi) ||
      !identical(nz$source$nodes, object$nodes) ||
      !identical(nz$source$centering, object$centering)) {
    aersn_abort("The supplied normalizer was built from different influence ",
                "contributions or a different variance-accumulation profile. ",
                "Build it from this fit (or pass a method name).")
  }
  check_method_grid(nz$method, object)
  invisible(nz)
}

check_method_grid <- function(method, object) {
  if (method_spec(method)$needs_calendar && object$centering != "calendar") {
    aersn_abort("Method '", method, "' requires calendar-time centering; ",
                "a supplied variance-accumulation profile is not supported.")
  }
  invisible(object)
}

#' @export
print.aersn_normalizer <- function(x, digits = 4L, ...) {
  cat(x$label, " normalizer (aersn_normalizer)\n", sep = "")
  cat("  method: ", x$method, "   family: ", x$family,
      "   statistic scale: ",
      if (x$scale == "gauge") "gauge (degree 1 in the estimation error)" else
        "Wald (degree 2 in the estimation error)", "\n", sep = "")
  cat("  n = ", x$n, "   q = ", x$q, "   parameters: ",
      paste(x$names, collapse = ", "), "\n", sep = "")
  if (length(x$tuning)) {
    cat("  tuning actually used:\n")
    for (nm in names(x$tuning)) {
      cat("    ", nm, " = ", format_tuning_value(x$tuning[[nm]]), "\n", sep = "")
    }
  } else {
    cat("  tuning: none (tuning-free)\n")
  }
  if (!is.null(x$matrix)) {
    cat("  normalizer matrix (", nrow(x$matrix), " x ", ncol(x$matrix),
        "):\n", sep = "")
    print(signif(x$matrix, digits))
    ev <- eigen((x$matrix + t(x$matrix)) / 2, symmetric = TRUE,
                only.values = TRUE)$values
    cat("    eigenvalues: ", paste(signif(ev, digits), collapse = ", "),
        "   condition ", fmt_num(max(ev) / min(ev), digits), "\n", sep = "")
  }
  if (!is.null(x$hull)) {
    cat("  increment hull: numerical rank ", x$hull$rank, " of ", x$q,
        ", condition ", fmt_num(x$hull$condition, digits), "\n", sep = "")
  }
  cat("  reference family: ", x$reference_family, "\n", sep = "")
  if (length(x$notes)) for (nt in x$notes) cat("  note: ", nt, "\n", sep = "")
  invisible(x)
}

format_tuning_value <- function(v) {
  if (is.null(v)) return("NULL")
  if (is.character(v) || is.logical(v)) return(paste(v, collapse = ", "))
  if (is.numeric(v) && all(is.finite(v)) && all(v == round(v)) &&
      all(abs(v) < 1e9)) {
    return(paste(format(as.integer(v)), collapse = ", "))
  }
  paste(vapply(v, function(z) format(signif(z, 6L)), character(1L)),
        collapse = ", ")
}

## ---------------------------------------------------------------------------
## The two operations every method must provide.
##
##   normalizer_statistic(nz, Z): the statistic at one or more values of
##     z = sqrt(n) (theta_hat - v), supplied as rows of Z.
##
##   normalizer_radius(nz, A, cv): for each row a of A, the half-width
##     sup{ a'(theta_hat - v) : T(v) <= cv }, already divided by sqrt(n).
## ---------------------------------------------------------------------------

normalizer_statistic <- function(nz, Z) UseMethod("normalizer_statistic")
normalizer_radius <- function(nz, A, cv) UseMethod("normalizer_radius")

#' @export
normalizer_statistic.aersn_normalizer <- function(nz, Z) {
  ## Default: quadratic family.
  Z <- as_direction_matrix(Z, nz$q)
  V <- nz$matrix
  ch <- quadratic_chol(V, nz)
  apply(Z, 1L, function(z) sum(backsolve(ch, z, transpose = TRUE)^2))
}

#' @export
normalizer_radius.aersn_normalizer <- function(nz, A, cv) {
  A <- as_direction_matrix(A, nz$q)
  V <- nz$matrix
  quadratic_chol(V, nz)
  quad <- rowSums((A %*% V) * A)
  if (any(quad < 0)) {
    aersn_abort("The normalizer matrix gave a negative quadratic form for a ",
                "contrast; it is not numerically positive semidefinite.")
  }
  sqrt(cv * quad) / sqrt(nz$n)
}

#' @export
normalizer_statistic.aersn_normalizer_hull <- function(nz, Z) {
  aersn_gauge(nz$hull, Z)
}

#' @export
normalizer_radius.aersn_normalizer_hull <- function(nz, A, cv) {
  cv / sqrt(nz$n) * aersn_support(nz$hull, A)
}

## Cholesky factor of a quadratic normalizer, with an informative failure.
quadratic_chol <- function(V, nz) {
  V <- (V + t(V)) / 2
  d <- sqrt(pmax(diag(V), 0))
  if (any(!is.finite(V)) || any(d <= 0)) {
    aersn_abort("The ", nz$label, " normalizer matrix is not positive definite.")
  }
  # Assess dependence after removing coordinate units, not on raw eigenvalues.
  C <- sweep(sweep(V, 1L, d, "/"), 2L, d, "/")
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(ev)) {
    aersn_abort("The ", nz$label, " normalizer matrix is numerically ",
                "singular after coordinate scaling. No statistic or joint ",
                "region is reported. No ridge term or pseudo-inverse is applied.")
  }
  sweep(chol(C), 2L, d, "*")
}

## Directions of the confidence region's support, used for plotting.
normalizer_boundary_2d <- function(nz, coords, n_angles) {
  phi <- seq(0, 2 * pi, length.out = n_angles + 1L)[-(n_angles + 1L)]
  D <- matrix(0, n_angles, nz$q)
  D[, coords[1L]] <- cos(phi)
  D[, coords[2L]] <- sin(phi)
  D
}

#' Hull normalizer
#'
#' Wraps the increment hull of an [aersn()] fit as an
#' [aersn_normalizer()] object, so that the affine-equivariant adjusted-range
#' method can be handled by the same inference and comparison functions as the
#' other methods.  The hull itself is computed by [aersn_hull()] when the fit
#' is created; this function does not recompute it.
#'
#' @param object An [aersn()] fit.
#' @return An `"aersn_normalizer"` object with `family = "geometric"`.
#' @seealso [aersn_normalizer()], [aersn_hull()].
#' @examples
#' set.seed(1)
#' fit <- aersn_mean(matrix(rnorm(200), 100, 2))
#' aersn_hull_normalizer(fit)
#' @export
aersn_hull_normalizer <- function(object) {
  check_aersn(object)
  new_normalizer("hull", q = object$q, n = object$n, names = object$names,
                 hull = object$hull,
                 diagnostics = list(rank = object$hull$rank,
                                    condition = object$hull$condition,
                                    min_projected_range =
                                      object$hull$min_projected_range),
                 centering = object$centering, nodes = object$nodes, source = object)
}
