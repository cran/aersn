#' Variance-accumulation profile on the observation grid
#'
#' Validates, or constructs from a function, the nodes
#' \eqn{\tau(k/n)}, \eqn{k = 0, \ldots, n}, of a variance-accumulation
#' profile used for profile centering of the influence path.  The profile
#' must be continuous and nondecreasing from \eqn{\tau(0) = 0} to
#' \eqn{\tau(1) = 1}; the package only checks these properties on the grid.
#'
#' @details
#' Under calendar-time centering the nodes are \eqn{k/n}.  The manuscript
#' (Section 4, "Nonuniform Variance Accumulation") replaces \eqn{k/n} by
#' \eqn{\tau(k/n)} when the cumulative long-run covariance satisfies the
#' fixed-shape condition \eqn{C(r) = \tau(r)\,C(1)}.  Two additional
#' conditions are required for a *feasible* profile-centered path
#' (manuscript Assumption "Feasible scalar time change"):
#' * the estimated profile belongs to
#'   \eqn{\{(t_0, \ldots, t_n): 0 = t_0 \le \cdots \le t_n = 1\}} and is
#'   uniformly consistent for \eqn{\tau}; and
#' * the increment hull of the feasible profile-centered path converges in
#'   Hausdorff distance to the hull of the infeasible one.
#'
#' Consistency of the profile alone does not deliver the second condition:
#' for the sample mean, the fitted contributions sum to zero, so
#' substituting any profile leaves the calendar-centered path unchanged,
#' and the resulting path differs from the required one by
#' \eqn{\{\tau(r) - r\}S_n(1)} (Supplement, Section "Cumulative information
#' and profile-based reindexing").  The manuscript establishes feasibility
#' for conditional-likelihood scores whose cumulative score variance and
#' expected negative Hessian accumulate the same fraction of their
#' full-sample matrices; see [aersn_mle()] and [aersn_opg_profile()].  The
#' package does not estimate a valid common profile for arbitrary dependent
#' data.
#'
#' @param profile Either `NULL` (calendar time, nodes `k/n`), a function of
#'   the sample fraction `r` in `[0, 1]`, or a numeric vector of length
#'   `n + 1` giving the nodes \eqn{\tau(k/n)} for \eqn{k = 0, \ldots, n}.
#' @param n Number of observations (the number of grid intervals).
#' @param tol Tolerance for the endpoint and monotonicity checks, between
#'   zero and `1e-6`.
#'
#' @return A numeric vector of length `n + 1` with attribute
#'   `"centering"` equal to `"calendar"` or `"profile"`.
#' @examples
#' aersn_profile(NULL, 4)
#' aersn_profile(function(r) r^3, 4)
#' aersn_profile(c(0, 0.1, 0.3, 0.7, 1), 4)
#' @export
aersn_profile <- function(profile, n, tol = 1e-10) {
  check_scalar_number(n, "n", lower = 1, integer = TRUE)
  check_scalar_number(tol, "tol", lower = 0, upper = 1e-6)
  n <- as.integer(n)
  grid <- seq.int(0L, n) / n
  if (is.null(profile)) {
    return(structure(grid, centering = "calendar"))
  }
  if (is.function(profile)) {
    nodes <- as.numeric(profile(grid))
    if (length(nodes) != n + 1L) {
      aersn_abort("`profile` as a function must return one value per grid ",
                  "point; got ", length(nodes), " values for ", n + 1L,
                  " points.")
    }
  } else if (is.numeric(profile)) {
    nodes <- as.numeric(profile)
    if (length(nodes) != n + 1L) {
      aersn_abort("`profile` must have length n + 1 = ", n + 1L,
                  " (the nodes tau(k/n) for k = 0, ..., n); got length ",
                  length(nodes), ".")
    }
  } else {
    aersn_abort("`profile` must be NULL, a function, or a numeric vector.")
  }
  if (any(!is.finite(nodes))) {
    aersn_abort("`profile` contains non-finite values.")
  }
  if (abs(nodes[1L]) > tol || abs(nodes[n + 1L] - 1) > tol) {
    aersn_abort("A variance-accumulation profile must satisfy tau(0) = 0 ",
                "and tau(1) = 1; got tau(0) = ", nodes[1L], " and tau(1) = ",
                nodes[n + 1L], ".")
  }
  d <- diff(nodes)
  if (any(d < -tol)) {
    k <- which(d < -tol)[1L]
    aersn_abort("A variance-accumulation profile must be nondecreasing; ",
                "it decreases between nodes ", k - 1L, " and ", k,
                " (from ", nodes[k], " to ", nodes[k + 1L], ").")
  }
  nodes <- cummax(pmin(pmax(nodes, 0), 1))
  nodes[1L] <- 0
  nodes[n + 1L] <- 1
  if (isTRUE(all.equal(nodes, grid, tolerance = 1e-12))) {
    return(structure(grid, centering = "calendar"))
  }
  structure(nodes, centering = "profile")
}

#' Centered influence path
#'
#' Computes the centered partial-sum path of estimated influence
#' contributions on the observation grid, including the initial and final
#' zero points.  This is the object whose increment hull normalizes the
#' estimation error.
#'
#' @details
#' With contributions \eqn{\hat\psi_{n,t} \in R^q} in the rows of `psi`
#' and profile nodes \eqn{\tau_k = \tau(k/n)}, the path is
#' \deqn{\hat G_n(k/n) = n^{-1/2}\Big\{\sum_{t=1}^k \hat\psi_{n,t}
#'   - \tau_k \sum_{t=1}^n \hat\psi_{n,t}\Big\}, \qquad k = 0, \ldots, n,}
#' manuscript equations (feasible-path) and (feasible-clock-path).  Under
#' calendar-time centering \eqn{\tau_k = k/n}.  The empty sum at \eqn{k = 0}
#' is zero, so \eqn{\hat G_n(0) = \hat G_n(1) = 0}.  The scaling by
#' \eqn{n^{-1/2}} uses `n = nrow(psi)`, which must be the number of
#' contributions in the asymptotically linear representation.
#'
#' Centering subtracts \eqn{\tau_k} times the full-sample sum, not the
#' sample mean of the contributions at each date; the two coincide only
#' under calendar-time centering.  If the contributions sum to zero (fitted
#' scores or least-squares moments), the second term vanishes and the
#' profile does not change the path; it still changes the reference grid
#' used by [aersn_reference()].
#'
#' @param psi Numeric vector (scalar parameter) or `n x q` numeric matrix of
#'   influence contributions with observations in rows and parameters in
#'   columns, in time order.
#' @param profile Variance-accumulation profile; see [aersn_profile()].
#'   `NULL` gives calendar-time centering.
#' @param names Optional parameter names (length `q`).
#'
#' @return An object of class `"aersn_path"`: a list with the
#'   `(n + 1) x q` matrix `G` of path values, `n`, `q`, the profile `nodes`,
#'   the `centering` label, and `total`, the full-sample sum of the
#'   contributions.
#' @seealso [aersn_hull()], [aersn()].
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(60), 30, 2)
#' path <- aersn_path(sweep(Y, 2, colMeans(Y)))
#' path$G[c(1, 31), ]   # zero at both ends
#' plot(path)
#' @export
aersn_path <- function(psi, profile = NULL, names = NULL) {
  psi <- as_contribution_matrix(psi, "psi")
  n <- nrow(psi)
  q <- ncol(psi)
  if (n < 2L) {
    aersn_abort("At least two observations are required; got n = ", n, ".")
  }
  names <- parameter_names(names, q, psi = psi)
  nodes <- aersn_profile(profile, n)
  centering <- attr(nodes, "centering")
  C <- apply(psi, 2L, cumsum)
  if (q == 1L) C <- matrix(C, ncol = 1L)
  total <- C[n, ]
  G <- rbind(rep(0, q), C - tcrossprod(as.numeric(nodes)[-1L], total)) /
    sqrt(n)
  G[n + 1L, ] <- 0   # exact: tau(1) = 1 removes the full-sample sum
  colnames(G) <- names
  rownames(G) <- NULL
  structure(list(G = G, n = n, q = q, nodes = as.numeric(nodes),
                 centering = centering, total = total, names = names),
            class = "aersn_path")
}

#' @export
print.aersn_path <- function(x, ...) {
  cat("Centered influence path (aersn_path)\n")
  cat("  n =", x$n, "observations,", x$q,
      if (x$q == 1L) "parameter" else "parameters",
      paste0("(", paste(x$names, collapse = ", "), ")"), "\n")
  cat("  centering:", if (x$centering == "calendar")
    "calendar time, tau(r) = r" else "supplied variance-accumulation profile",
    "\n")
  cat("  path values: (n + 1) x q matrix `G`, G[1, ] = G[n + 1, ] = 0\n")
  cat("  coordinate adjusted ranges:",
      paste(fmt_num(column_ranges(x$G)), collapse = ", "), "\n")
  invisible(x)
}

#' Feasible variance-accumulation profile from score outer products
#'
#' Estimates the variance-accumulation profile \eqn{\tau} from fitted
#' conditional-likelihood scores by the normalized trace of cumulative
#' outer products, manuscript equation (opg-trace-clock):
#' \deqn{\hat\tau_n(k/n) = p^{-1}\,\mathrm{tr}\big(\hat\Sigma_{s,n}^{-1/2}
#'   \hat\Sigma_{s,k}\hat\Sigma_{s,n}^{-1/2}\big), \qquad
#'   \hat\Sigma_{s,k} = n^{-1}\sum_{t \le k} s_t s_t^\top .}
#' The estimate is nondecreasing, lies in \eqn{[0, 1]}, and has endpoints
#' zero and one whenever \eqn{\hat\Sigma_{s,n}} is positive definite.
#'
#' @details
#' The estimator is justified (Supplement, Proposition "Feasible
#' variance-accumulation profile from score outer products") for
#' martingale-difference conditional scores whose cumulative conditional
#' variance and cumulative expected negative Hessian accumulate the same
#' fraction \eqn{\tau(r)} of their full-sample matrices.  For serially
#' correlated influence contributions, lag-zero outer products do not
#' estimate cumulative long-run covariance, and this profile is not valid.
#' The returned proportionality measure
#' \eqn{\mathcal E_n = \max_k p^{-1/2}\|\hat\Xi_{s,k} - \hat\tau_n(k/n)I_p\|_F}
#' (manuscript equation (clock-proportionality-measure)) is descriptive:
#' large values indicate that the accumulation is not proportional across
#' score components, so a common scalar profile is doubtful.
#'
#' @param scores `n x p` matrix of fitted conditional scores
#'   \eqn{s_t(\hat\beta_n)} in time order (all model parameters, including
#'   nuisance parameters).
#'
#' @return A numeric vector of length `n + 1` (the nodes, with attribute
#'   `"centering" = "profile"`), with attributes `"proportionality"`
#'   (\eqn{\mathcal E_n}) and `"score_sum"` (the scaled full-sample score
#'   \eqn{n^{-1/2}\sum_t s_t}, which should be close to zero for fitted
#'   scores).
#' @seealso [aersn_mle()], [aersn_profile()].
#' @examples
#' set.seed(2)
#' s <- matrix(rnorm(400), 200, 2)      # synthetic scores
#' s <- sweep(s, 2, colMeans(s))         # fitted scores sum to zero
#' tau <- aersn_opg_profile(s)
#' range(diff(tau)); attr(tau, "proportionality")
#' @export
aersn_opg_profile <- function(scores) {
  s <- as_contribution_matrix(scores, "scores")
  n <- nrow(s)
  p <- ncol(s)
  ## Diagonal scaling changes the whitening only by an orthogonal matrix;
  ## both the trace profile and Frobenius proportionality measure are
  ## unchanged. It avoids treating different score units as singularity.
  units <- apply(abs(s), 2L, max)
  if (any(units == 0)) aersn_abort("The score covariance is not positive definite.")
  ss <- sweep(s, 2L, units, `/`)
  Sigma_n <- crossprod(ss) / n
  W <- symmetric_inverse_sqrt(Sigma_n, "crossprod(scores)/n")
  s_white <- ss %*% W
  w <- rowSums(s_white^2) / (n * p)     # increments of tau
  tau <- c(0, cumsum(w))
  tau[n + 1L] <- 1
  ## Proportionality measure: max_k p^{-1/2} || Xi_k - tau_k I ||_F
  Xi <- matrix(0, p, p)
  prop <- 0
  for (k in seq_len(n)) {
    Xi <- Xi + tcrossprod(s_white[k, ]) / n
    dev <- Xi - diag(tau[k + 1L], p)
    prop <- max(prop, sqrt(sum(dev^2)) / sqrt(p))
  }
  structure(cummax(pmin(pmax(tau, 0), 1)), centering = "profile",
            proportionality = prop,
            score_sum = colSums(s) / sqrt(n))
}
