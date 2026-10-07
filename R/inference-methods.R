## Method dispatch shared by aersn_test(), confint(), aersn_contrast(),
## aersn_region() and aersn_compare().
##
## Inference on a linear target is always carried out by rebuilding the
## method on that target's influence contributions, rather than by projecting
## the full-dimensional normalizer. The two agree for the increment hull and
## for every method whose normalizer is a linear functional of outer products,
## and they differ, deliberately, where the method itself is not: the LDL
## factorization is recomputed for the target, and a data-dependent HAC
## bandwidth is re-selected from the target's contributions.

## Collect method tuning from `...` together with the arguments that are
## given their own formal.
##
## `nu` needs a formal of its own: it is a prefix of `null`, so without an
## exact formal R's partial matching would silently send `nu = 16` to the null
## value. Exact matching on a formal takes precedence over partial matching,
## which removes the trap. Any other tuning name that becomes a prefix of a
## formal must be given the same treatment; the guard below detects the case.
method_tuning <- function(dots, nu = NULL, call = NULL) {
  if (!is.null(nu)) {
    if (!is.null(dots$nu)) {
      aersn_abort("`nu` was supplied more than once.")
    }
    dots$nu <- nu
  }
  if (!is.null(call)) {
    supplied <- names(call)
    supplied <- supplied[!is.na(supplied) & nzchar(supplied)]
    absorbed <- setdiff(intersect(supplied, normalizer_tuning_names()),
                        names(dots))
    if (length(absorbed)) {
      aersn_abort("Tuning argument(s) ", paste(absorbed, collapse = ", "),
                  " were absorbed by partial matching of another argument ",
                  "name and did not reach the method. Build the normalizer ",
                  "explicitly with aersn_normalizer() and pass it as ",
                  "`method`.")
    }
  }
  dots
}

## Resolve `method` (a string with tuning in `tuning`, or a prebuilt object).
resolve_method <- function(object, method, tuning = list()) {
  if (inherits(method, "aersn_normalizer")) {
    if (length(tuning)) {
      aersn_abort("Tuning arguments (", paste(names(tuning), collapse = ", "),
                  ") cannot be combined with a prebuilt normalizer. Either ",
                  "pass method = \"", method$method, "\" with the tuning ",
                  "arguments, or build the normalizer with them.")
    }
    check_normalizer_fit(method, object)
    return(method)
  }
  do.call(aersn_normalizer, c(list(object, method), tuning))
}

## Build the fit for a linear target, or return the fit itself for the
## identity transformation.
target_fit <- function(object, A, labels) {
  if (nrow(A) == object$q &&
      isTRUE(all.equal(unname(A), unname(diag(object$q)), tolerance = 1e-12))) {
    return(object)
  }
  aersn_target(object, A, names = labels)
}

## Everything the inference functions need for one set of linear contrasts.
##
##   A       m x q contrast matrix in the original parameter coordinates
##   type    "simultaneous" (project the q-dimensional region),
##           "joint"        (m-dimensional region for these contrasts),
##           "marginal"     (a separate one-dimensional problem per contrast)
method_intervals <- function(object, A, labels, type, method, tuning, level,
                             reference, draws, seed, cores) {
  m <- nrow(A)
  est <- as.numeric(A %*% object$estimate)
  if (type == "simultaneous") {
    nz <- resolve_method(object, method, tuning)
    require_method_full_rank(object, nz, "simultaneous inference")
    ref <- reference_for_normalizer(nz, object$q, reference, object$n,
                                    object$nodes,
                                    object$centering == "calendar",
                                    draws, seed, cores)
    cv <- aersn_critical_value(ref, level)
    half <- normalizer_radius(nz, A, cv)
    return(list(estimate = est, half = half, cv = cv, dim = object$q,
                nz = nz, ref = ref))
  }
  if (type == "joint") {
    if (qr(A)$rank < m) {
      aersn_abort("`type = \"joint\"` requires contrasts of full row rank.")
    }
    tf <- target_fit(object, A, labels)
    nz <- resolve_method(tf, method, tuning)
    require_method_full_rank(tf, nz, "joint inference on the selected target")
    ref <- reference_for_normalizer(nz, m, reference, tf$n, tf$nodes,
                                    tf$centering == "calendar",
                                    draws, seed, cores)
    cv <- aersn_critical_value(ref, level)
    half <- normalizer_radius(nz, diag(m), cv)
    return(list(estimate = est, half = half, cv = cv, dim = m,
                nz = nz, ref = ref))
  }
  ## marginal: one one-dimensional problem per contrast, sharing one reference
  nz1 <- NULL
  ref <- NULL
  cv <- NA_real_
  half <- numeric(m)
  normalizers <- stats::setNames(vector("list", m), labels)
  for (i in seq_len(m)) {
    tf <- target_fit(object, A[i, , drop = FALSE], labels[i])
    nzi <- resolve_method(tf, method, tuning)
    normalizers[[i]] <- nzi
    if (i == 1L) {
      nz1 <- nzi
      ref <- reference_for_normalizer(nzi, 1L, reference, tf$n, tf$nodes,
                                      tf$centering == "calendar",
                                      draws, seed, cores)
      cv <- aersn_critical_value(ref, level)
    }
    half[i] <- normalizer_radius(nzi, matrix(1, 1L, 1L), cv)
  }
  list(estimate = est, half = half, cv = cv, dim = 1L, nz = nz1, ref = ref,
       tuning_by_contrast = lapply(normalizers, `[[`, "tuning"))
}

## Rank requirement: the hull needs a full-dimensional increment hull, the
## quadratic methods need a nonsingular normalizer matrix. The matrix is
## checked where it is used, which reports the eigenvalues.
require_method_full_rank <- function(object, nz, what) {
  if (nz$family == "geometric") require_full_rank(object, what)
  else quadratic_chol(nz$matrix, nz)
  invisible(nz)
}

## Selection matrix for coordinates.
selection_matrix <- function(idx, q) {
  A <- diag(q)[idx, , drop = FALSE]
  A
}

## Label for the printed method line.
method_line <- function(nz) {
  tun <- nz$tuning
  if (!length(tun)) return(nz$label)
  shown <- vapply(names(tun), function(nm) {
    paste0(nm, " = ", format_tuning_value(tun[[nm]]))
  }, character(1L))
  paste0(nz$label, " (", paste(shown, collapse = ", "), ")")
}

#' Compare inference methods on one estimate
#'
#' Applies several methods to the same estimate, the same influence
#' contributions, the same null value or linear restriction and the same
#' confidence level, and collects the results in one table.  Each row reports
#' the statistic, its critical value, a p-value, the reference law and the
#' tuning values actually used, so that differences between methods can be
#' attributed to the normalizer rather than to a difference in setup.
#'
#' @details
#' Statistics are not comparable across rows as numbers.  The increment-hull
#' statistic is a gauge, homogeneous of degree one in the estimation error,
#' while the other five are Wald statistics, homogeneous of degree two, and
#' each is referred to its own law.  The comparable quantities are the
#' decisions, the p-values and the interval widths, which the `reject`,
#' `p_value` and `half_width` columns report.  The `half_width` column is the
#' half-width of the confidence interval for the first contrast, or for the
#' first coordinate when no contrast is given, at the stated level.
#'
#' Methods that fail on the data, for example because a normalizer matrix is
#' numerically singular, are reported with `NA` statistics and the reason in
#' the `status` column rather than being dropped from the table.
#'
#' @param object An [aersn()] fit.
#' @param null Null value for the point-null test; defaults to zeros.
#' @param contrast Optional `m` by `q` matrix of full row rank.  When
#'   supplied the comparison is of the restriction `contrast %*% theta = null`
#'   and `null` has length `m`.
#' @param level Confidence level.
#' @param methods Method identifiers; see [aersn_methods()].  The default is
#'   all six, reported in the order increment hull, LDL componentwise,
#'   quadratic self-normalization, HAC, Bartlett fixed-b, EWC.
#' @param settings A named list of per-method tuning arguments, for example
#'   `list(hac = list(kernel = "Parzen", bandwidth = "andrews"),
#'   fixedb = list(b = 0.5), ewc = list(nu = 16))`.
#' @param draws,seed,cores Settings for any simulated reference law.
#' @param verbose Report progress while reference laws are simulated.
#'
#' @return A data frame of class `"aersn_comparison"` with one row per method
#'   and columns `method`, `label`, `statistic`, `statistic_scale`,
#'   `critical_value`, `p_value`, `p_mcse`, `reject`, `half_width`,
#'   `reference_family`, `reference`, `tuning` and `status`.  The estimate,
#'   null value and level are kept as attributes.
#' @seealso [aersn_test()], [aersn_normalizer()].
#' @examples
#' set.seed(20)
#' n <- 300
#' e <- matrix(rnorm(2 * n), n, 2)
#' Y <- e
#' for (t in 2:n) Y[t, ] <- 0.4 * Y[t - 1, ] + e[t, ]
#' fit <- aersn_mean(Y, names = c("m1", "m2"))
#' aersn_compare(fit, null = c(0, 0), draws = 400, seed = 1)
#' @export
aersn_compare <- function(object, null = NULL, contrast = NULL, level = 0.95,
                          methods = aersn_methods(), settings = list(),
                          draws = NULL, seed = 1L, cores = 1L,
                          verbose = FALSE) {
  check_aersn(object)
  check_level(level)
  check_flag(verbose, "verbose")
  if (!is.character(methods) || !length(methods)) {
    aersn_abort("`methods` must be a character vector of method identifiers.")
  }
  if (anyDuplicated(methods)) {
    aersn_abort("`methods` must not repeat a method.")
  }
  unknown <- setdiff(methods, aersn_methods())
  if (length(unknown)) {
    aersn_abort("Unknown method(s): ", paste(unknown, collapse = ", "),
                ". Supported methods are ",
                paste(aersn_methods(), collapse = ", "), ".")
  }
  ## Keep the documented reporting order whatever order the user supplies.
  methods <- aersn_methods()[aersn_methods() %in% methods]
  if (!is.list(settings) ||
      (length(settings) && (is.null(names(settings)) || anyNA(names(settings)) ||
                            any(!nzchar(names(settings))) || anyDuplicated(names(settings))))) {
    aersn_abort("`settings` must be a named list of per-method tuning lists.")
  }
  bad <- setdiff(names(settings), aersn_methods())
  if (length(bad)) {
    aersn_abort("`settings` has entries for unknown method(s): ",
                paste(bad, collapse = ", "), ".")
  }
  work <- object
  if (!is.null(contrast)) {
    A <- as_direction_matrix(contrast, object$q)
    labels <- if (!is.null(rownames(A))) rownames(A) else
      paste0("contrast", seq_len(nrow(A)))
    if (qr(A)$rank < nrow(A)) {
      aersn_abort("`contrast` must have full row rank.")
    }
    work <- aersn_target(object, A, names = labels)
  }
  q <- work$q
  if (is.null(null)) null <- rep(0, q)
  null <- as.numeric(null)
  if (length(null) != q || any(!is.finite(null))) {
    aersn_abort("`null` must be a finite vector of length ", q, ".")
  }
  names(null) <- work$names
  rows <- lapply(methods, function(mth) {
    tuning <- if (is.null(settings[[mth]])) list() else settings[[mth]]
    if (!is.list(tuning)) {
      aersn_abort("`settings$", mth, "` must be a list of tuning arguments.")
    }
    if (verbose) message("Comparing method: ", mth)
    out <- tryCatch({
      nz <- resolve_method(work, mth, tuning)
      require_method_full_rank(work, nz, "the comparison")
      z <- sqrt(work$n) * (work$estimate - null)
      stat <- normalizer_statistic(nz, matrix(z, nrow = 1L))
      ref <- reference_for_normalizer(nz, q, NULL, work$n, work$nodes,
                                      work$centering == "calendar",
                                      draws, seed, cores, verbose = verbose)
      cv <- aersn_critical_value(ref, level)
      pv <- aersn_pvalue(ref, stat)
      hw <- normalizer_radius(nz, matrix(diag(q)[1L, ], nrow = 1L), cv)
      data.frame(method = mth, label = nz$label,
                 statistic = stat,
                 statistic_scale = nz$scale,
                 critical_value = cv,
                 p_value = pv$p.value,
                 p_mcse = if (is.null(pv$mcse)) NA_real_ else pv$mcse,
                 reject = stat > cv,
                 half_width = hw,
                 reference_family = nz$reference_family,
                 reference = reference_summary(ref),
                 tuning = tuning_string(nz$tuning),
                 status = "ok",
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(method = mth, label = method_spec(mth)$label,
                 statistic = NA_real_, statistic_scale = method_spec(mth)$scale,
                 critical_value = NA_real_, p_value = NA_real_,
                 p_mcse = NA_real_, reject = NA, half_width = NA_real_,
                 reference_family = method_spec(mth)$reference,
                 reference = NA_character_, tuning = tuning_string(tuning),
                 status = conditionMessage(e), stringsAsFactors = FALSE)
    })
    out
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  structure(tab, class = c("aersn_comparison", "data.frame"),
            estimate = work$estimate, null.value = null, level = level,
            n = work$n, q = q,
            contrast = if (is.null(contrast)) NULL else contrast)
}

tuning_string <- function(tuning) {
  if (!length(tuning)) return("none")
  paste(vapply(names(tuning), function(nm) {
    paste0(nm, "=", format_tuning_value(tuning[[nm]]))
  }, character(1L)), collapse = "; ")
}

#' @export
`[.aersn_comparison` <- function(x, ...) {
  drop_result_class(NextMethod(), x, "aersn_comparison",
                    c("estimate", "null.value", "level", "n", "q", "contrast"))
}

#' @export
print.aersn_comparison <- function(x, digits = 4L, ...) {
  if (is.null(attr(x, "n")) || is.null(attr(x, "estimate"))) {
    return(print(as.data.frame(unclass(x)), digits = digits))
  }
  cat("Comparison of inference methods on one estimate\n")
  cat("  n = ", attr(x, "n"), "   q = ", attr(x, "q"),
      "   level = ", attr(x, "level"), "\n", sep = "")
  cat("  null value: ",
      paste(names(attr(x, "null.value")), "=",
            fmt_num(attr(x, "null.value"), digits), collapse = ", "), "\n",
      sep = "")
  cat("  estimate:   ",
      paste(names(attr(x, "estimate")), "=",
            fmt_num(attr(x, "estimate"), digits), collapse = ", "), "\n",
      sep = "")
  show <- data.frame(
    method = x$method,
    statistic = ifelse(is.na(x$statistic), "-", fmt_num(x$statistic, digits)),
    scale = x$statistic_scale,
    critical = ifelse(is.na(x$critical_value), "-",
                      fmt_num(x$critical_value, digits)),
    p = ifelse(is.na(x$p_value), "-", fmt_num(x$p_value, 3)),
    reject = ifelse(is.na(x$reject), "-", ifelse(x$reject, "yes", "no")),
    half_width = ifelse(is.na(x$half_width), "-",
                        fmt_num(x$half_width, digits)),
    stringsAsFactors = FALSE)
  print(show, row.names = FALSE, right = FALSE)
  cat("\nStatistics are on different scales and are not comparable as numbers:",
      "\n  the hull statistic is a gauge, the others are Wald statistics,",
      "\n  and each is referred to its own law. The comparable columns are",
      "\n  p, reject and half_width (for the first coordinate).\n")
  cat("\nTuning actually used:\n")
  for (i in seq_len(nrow(x))) {
    cat("  ", format(x$method[i], width = 7), x$tuning[i], "\n", sep = "")
  }
  failed <- x$status != "ok"
  if (any(failed)) {
    cat("\nMethods that could not be evaluated:\n")
    for (i in which(failed)) {
      cat("  ", x$method[i], ": ", x$status[i], "\n", sep = "")
    }
  }
  invisible(x)
}
