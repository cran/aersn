#' Increment-hull self-normalizer
#'
#' Constructs the affine-equivariant normalizing set of the manuscript:
#' the convex hull of all increments of the centered influence path,
#' \deqn{K(\hat G_n) = \mathrm{conv}\{\hat G_n(r) - \hat G_n(s):
#'   0 \le r, s \le 1\},}
#' manuscript equation (increment-hull).  The set is compact, convex and
#' centrally symmetric.  Its support function in direction \eqn{u} is the
#' adjusted range of the projected path,
#' \eqn{h_K(u) = \max_j u^\top \hat G_{n,j} - \min_j u^\top \hat G_{n,j}}
#' (equation (support-range)), and its gauge has the dual representation
#' \eqn{\gamma_K(x) = \sup_{u \ne 0} |u^\top x| / h_K(u)} (equation
#' (gauge-dual)).  The hull is never stored as a list of pairwise
#' increments; it is represented by the path itself, and the gauge is
#' evaluated by the linear program of equation (gauge-lp).
#'
#' @details
#' The object records diagnostics recommended in the Supplement (Section
#' "Multivariate critical values and affine transformations"): the
#' singular values of the `(n + 1) x q` path matrix, its numerical rank
#' after dividing each nonconstant coordinate by its adjusted range
#' (threshold `max(n + 1, q) * s_1 * eps` on this scaled path), the
#' condition number `s_1 / s_q`, and the minimum, maximum and spread of
#' the projected adjusted range over the coordinate axes and a fixed set of
#' `n_directions` prespecified unit directions.  If the numerical rank is
#' below `q` the hull is not full dimensional and joint statistics are not
#' evaluated (Algorithm 1, step 2). Scaling prevents units alone from
#' causing numerical rank deficiency; unscaled diagnostics remain available
#' as `singular_values`, `raw_rank` and `condition`. The scanned minimum range is an upper
#' bound on the all-direction minimum.
#'
#' @param path An [aersn_path()] object.
#' @param n_directions Number of prespecified unit directions used for the
#'   projected-range diagnostics (in addition to the coordinate axes).
#'
#' @return An object of class `"aersn_hull"` with elements `path`, `G`,
#'   `n`, `q`, `singular_values`, `rank`, `full_rank`, `condition`,
#'   `coordinate_ranges`, `min_projected_range`, `max_projected_range`,
#'   `projected_range_spread`, and `n_directions`.
#' @seealso [aersn_gauge()], [aersn_support()], [aersn_path()].
#' @references Hong, Lin, Linton, Newey and Sun (2026), Sections 2.2-2.3
#'   and Supplement Section "Multivariate critical values and affine
#'   transformations".
#' @examples
#' set.seed(3)
#' Y <- matrix(rnorm(300), 100, 3)
#' hull <- aersn_hull(aersn_path(sweep(Y, 2, colMeans(Y))))
#' hull
#' aersn_support(hull, c(1, 0, 0))             # adjusted range of column 1
#' aersn_gauge(hull, c(0.5, -0.2, 0.1))
#' @export
aersn_hull <- function(path, n_directions = 1000L) {
  if (!inherits(path, "aersn_path")) {
    aersn_abort("`path` must be an aersn_path object; see aersn_path().")
  }
  check_scalar_number(n_directions, "n_directions", lower = 0,
                      integer = TRUE)
  G <- path$G
  n <- path$n
  q <- path$q
  if (n < q + 1L) {
    aersn_abort("At least q + 1 = ", q + 1L, " observations are needed for ",
                "a full-dimensional increment hull; got n = ", n, ".")
  }
  sv <- svd(G, nu = 0L, nv = 0L)$d
  raw_tol <- max(dim(G)) * sv[1L] * .Machine$double.eps
  coord <- column_ranges(G)
  scaled_G <- sweep(G, 2L, ifelse(coord > 0, coord, 1), `/`)
  scaled_sv <- svd(scaled_G, nu = 0L, nv = 0L)$d
  tol <- max(dim(G)) * scaled_sv[1L] * .Machine$double.eps
  rank <- sum(scaled_sv > tol)
  full_rank <- rank == q
  condition <- if (full_rank && sv[q] > 0) sv[1L] / sv[q] else Inf
  names(coord) <- path$names
  if (q >= 2L) {
    dirs <- rbind(diag(q), halton_directions(as.integer(n_directions), q))
    ## Projected ranges are computed in blocks of directions so that the
    ## memory use stays proportional to n for long series.
    proj_ranges <- numeric(nrow(dirs))
    block <- max(1L, floor(2e6 / nrow(G)))
    for (start in seq.int(1L, nrow(dirs), by = block)) {
      idx <- start:min(nrow(dirs), start + block - 1L)
      proj <- dirs[idx, , drop = FALSE] %*% t(G)
      proj_ranges[idx] <- apply(proj, 1L, function(x) max(x) - min(x))
    }
  } else {
    dirs <- matrix(1, 1L, 1L)
    proj_ranges <- coord
  }
  structure(list(path = path, G = G, n = n, q = q, names = path$names,
                 singular_values = sv, rank = rank, full_rank = full_rank,
                 scaled_singular_values = scaled_sv,
                 raw_rank = sum(sv > raw_tol), raw_rank_tolerance = raw_tol,
                 rank_tolerance = tol, condition = condition,
                 coordinate_ranges = coord,
                 min_projected_range = min(proj_ranges),
                 max_projected_range = max(proj_ranges),
                 projected_range_spread = max(proj_ranges) / min(proj_ranges),
                 n_directions = nrow(dirs)),
            class = "aersn_hull")
}

#' @export
print.aersn_hull <- function(x, ...) {
  cat("Increment-hull self-normalizer (aersn_hull)\n")
  cat("  n =", x$n, " q =", x$q, " centering:", x$path$centering, "\n")
  cat("  numerical rank:", x$rank, "of", x$q,
      if (x$full_rank) "(full dimensional)" else
        "(RANK DEFICIENT: joint statistics are not evaluated)", "\n")
  cat("  SVD condition number of the centered path:",
      fmt_num(x$condition), "\n")
  cat("  coordinate adjusted ranges:",
      paste(fmt_num(x$coordinate_ranges), collapse = ", "), "\n")
  cat("  projected range over", x$n_directions, "directions: min",
      fmt_num(x$min_projected_range), " max", fmt_num(x$max_projected_range),
      " spread", fmt_num(x$projected_range_spread), "\n")
  invisible(x)
}

#' Support function of an increment hull or confidence region
#'
#' For an [aersn_hull()] object, returns the support function
#' \eqn{h_K(u) = \sup_{z \in K} u^\top z}, which equals the adjusted range of
#' the projected path \eqn{u^\top \hat G_n} (manuscript equation
#' (support-range)).  For an [aersn_region()] object, returns the support
#' function of the confidence region
#' \eqn{u^\top\hat\theta_n + (c/\sqrt n)\,h_K(u)}, which is the upper end of
#' the simultaneous projection interval for the contrast \eqn{u^\top\theta}.
#'
#' @param x An `"aersn_hull"` or `"aersn_region"` object.
#' @param u A direction vector of length `q`, or a matrix with one direction
#'   per row.
#' @param ... Unused.
#' @return A numeric vector with one value per direction.
#' @export
aersn_support <- function(x, u, ...) UseMethod("aersn_support")

#' @rdname aersn_support
#' @export
aersn_support.aersn_hull <- function(x, u, ...) {
  U <- as_direction_matrix(u, x$q)
  proj <- U %*% t(x$G)
  apply(proj, 1L, function(v) max(v) - min(v))
}

as_direction_matrix <- function(u, q) {
  if (!is.numeric(u) || length(dim(u)) > 2L) {
    aersn_abort("`u` must be a numeric vector or matrix.")
  }
  if (is.null(dim(u))) {
    if (length(u) != q) {
      aersn_abort("`u` has length ", length(u), " but q = ", q, ".")
    }
    U <- matrix(as.numeric(u), nrow = 1L)
  } else {
    U <- as.matrix(u)
    if (ncol(U) != q) {
      aersn_abort("`u` has ", ncol(U), " columns but q = ", q, ".")
    }
    storage.mode(U) <- "double"
  }
  if (any(!is.finite(U))) aersn_abort("`u` must be finite.")
  U
}

#' Gauge of an increment hull, or the adjusted-range statistic
#'
#' The gauge \eqn{\gamma_K(x) = \inf\{\lambda > 0: x \in \lambda K\}} is
#' the factor by which the increment hull must be enlarged to contain
#' \eqn{x}.  It is evaluated by the linear program of manuscript equation
#' (gauge-lp), which uses only the \eqn{n + 1} path values,
#' \deqn{\gamma_K(x) = \max_{u, \ell}\ u^\top x \quad\text{subject to}\quad
#'   \ell \le u^\top \hat G_{n,j} \le \ell + 1,\ j = 0, \ldots, n.}
#' For `q = 1` the gauge is \eqn{|x|} divided by the adjusted range of the
#' path (exact scalar reduction, equation (scalar-reduction)).
#'
#' Applied to an [aersn()] object with a candidate parameter value `v`, the
#' function returns the test statistic
#' \eqn{T_n(v) = \gamma_{K(\hat G_n)}\{\sqrt n(\hat\theta_n - v)\}}.
#' Applied to an [aersn_region()] object it returns the gauge relative to
#' the region, so that values at most one indicate membership.
#'
#' @param x An `"aersn_hull"`, `"aersn"` or `"aersn_region"` object.
#' @param v For a hull, the point(s) at which to evaluate the gauge.  For a
#'   fit or region, the candidate parameter value(s).  A vector of length
#'   `q` or a matrix with one candidate per row.
#' @param details If `TRUE`, return a list with the value, the maximizing
#'   direction `u` of the dual program (normalized so that the projected
#'   range in that direction is one), and solver information; only for a
#'   single candidate.
#' @param ... Passed to methods.
#' @return A numeric vector of gauge values, or a list when
#'   `details = TRUE`.
#' @examples
#' set.seed(4)
#' Y <- matrix(rnorm(200), 100, 2)
#' fit <- aersn_mean(Y)
#' aersn_gauge(fit, c(0, 0))                    # statistic for H0: mean = 0
#' aersn_gauge(fit, c(0, 0), details = TRUE)$direction
#' @export
aersn_gauge <- function(x, v, ...) UseMethod("aersn_gauge")

#' @rdname aersn_gauge
#' @export
aersn_gauge.aersn_hull <- function(x, v, details = FALSE, ...) {
  check_flag(details, "details")
  if (!x$full_rank) {
    aersn_abort("The increment hull is not full dimensional (numerical rank ",
                x$rank, " < q = ", x$q, "); the joint statistic is not ",
                "evaluated. Remove redundant parameters or supply ",
                "contributions that span R^q.")
  }
  V <- as_direction_matrix(v, x$q)
  if (details) {
    if (nrow(V) != 1L) aersn_abort("`details = TRUE` requires one point.")
    return(gauge_lp(x$G, V[1L, ]))
  }
  vapply(seq_len(nrow(V)), function(i) gauge_lp(x$G, V[i, ])$value,
         numeric(1L))
}
