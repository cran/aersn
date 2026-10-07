## Componentwise adjusted ranges after lag-zero LDL partial prewhitening.

#' LDL factorization of a symmetric positive-definite matrix
#'
#' Factorizes \eqn{S = L D L^\top} with \eqn{L} unit lower triangular and
#' \eqn{D} diagonal with positive entries.  The factorization is computed from
#' the Cholesky factor as \eqn{L = C \mathrm{diag}(C)^{-1}} and
#' \eqn{D = \mathrm{diag}(C)^2}, where \eqn{C} is the lower Cholesky factor,
#' so that \eqn{C = L D^{1/2}}.
#'
#' @param S A symmetric positive-definite matrix.
#' @param what A label used in error messages.
#' @return A list with `L` (unit lower triangular), `D` (the diagonal of
#'   \eqn{D}, a vector), `chol` (the lower Cholesky factor \eqn{C}) and
#'   `reconstruction_error`, the maximum absolute deviation of
#'   \eqn{L D L^\top} from the symmetrized input.
#' @details
#' No ridge term, pseudo-inverse or nearest-positive-definite projection is
#' applied.  If the matrix is not numerically positive definite the function
#' stops and reports the smallest eigenvalue, because a silently modified
#' factorization would change the statistic without changing its name.
#' @examples
#' S <- crossprod(matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 4), 3, 3))
#' f <- aersn_ldl(S)
#' f$L
#' max(abs(f$L %*% diag(f$D) %*% t(f$L) - S))
#' max(abs(f$chol - f$L %*% diag(sqrt(f$D))))
#' @export
aersn_ldl <- function(S, what = "matrix") {
  S <- as.matrix(S)
  storage.mode(S) <- "double"
  if (!nrow(S) || nrow(S) != ncol(S)) {
    aersn_abort("`", what, "` must be square; it is ", nrow(S), " by ",
                ncol(S), ".")
  }
  if (any(!is.finite(S))) aersn_abort("`", what, "` must be finite.")
  Ssym <- (S + t(S)) / 2
  if (max(abs(S - Ssym)) > 1e-8 * max(1, max(abs(S)))) {
    aersn_abort("`", what, "` is not symmetric.")
  }
  dscale <- sqrt(pmax(diag(Ssym), 0))
  if (any(dscale <= 0)) {
    aersn_abort("`", what, "` is not positive definite.")
  }
  scaled <- sweep(sweep(Ssym, 1L, dscale, "/"), 2L, dscale, "/")
  ev <- eigen(scaled, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(ev)) {
    aersn_abort("`", what, "` is not numerically positive definite after ",
                "coordinate scaling. No ridge term is applied.")
  }
  C <- sweep(t(chol(scaled)), 1L, dscale, "*")
  d <- diag(C)
  L <- C %*% diag(1 / d, nrow = length(d))
  D <- d^2
  err <- max(abs(L %*% diag(D, nrow = length(D)) %*% t(L) - Ssym))
  list(L = L, D = D, chol = C, reconstruction_error = err)
}

#' Componentwise adjusted-range normalizer after LDL partial prewhitening
#'
#' Implements the componentwise adjusted-range self-normalizer of
#' Hong, Linton, McCabe, Sun and Wang (2024), as described in the
#' Supplement of the package's main reference.  The influence contributions
#' are transformed by the inverse of the unit lower triangular factor of
#' their sample lag-zero covariance, which removes contemporaneous cross
#' dependence, and each transformed coordinate is normalized by its own
#' adjusted range.  The joint statistic is the sum of the squared
#' componentwise ratios.
#'
#' @details
#' Write \eqn{\hat\Gamma_0} for the sample lag-zero covariance of the
#' influence contributions, factorized as
#' \eqn{\hat\Gamma_0 = \hat L \hat D \hat L^\top} with \eqn{\hat L} unit lower
#' triangular ([aersn_ldl()]).  With \eqn{z = \sqrt n(\hat\theta_n - v)},
#' transformed error \eqn{\tilde z = \hat L^{-1} z} and transformed path
#' \eqn{\tilde G_n = \hat L^{-1}\hat G_n}, the statistic is
#' \deqn{T^{\mathrm{comp}}_n(v) = \sum_{j=1}^q
#'   \left(\frac{|\tilde z_j|}{\max_k \tilde G_{n,k,j}
#'     - \min_k \tilde G_{n,k,j}}\right)^2 ,}
#' the range being taken over the observation grid including the origin.
#' This is the Wald form \eqn{z^\top V^{-1} z} with
#' \eqn{V = \hat L \,\mathrm{diag}(R_j^2)\, \hat L^\top}, where \eqn{R_j} is
#' the adjusted range of transformed coordinate \eqn{j}, so the confidence
#' region is an ellipsoid in the original parameter coordinates and linear
#' contrasts are obtained from \eqn{V} in the usual way.
#'
#' **Cholesky and LDL give the same statistic.** Writing the lower Cholesky
#' factor as \eqn{\hat C = \hat L \hat D^{1/2}}, the transformed error becomes
#' \eqn{\hat D^{-1/2}\tilde z} and transformed coordinate \eqn{j} of the path
#' is divided by the same \eqn{\hat D_{jj}^{1/2}}, so each ratio is unchanged
#' and the diagonal scale cancels exactly.  The function returns both factors
#' and reports the difference between the two statistics as a diagnostic.
#'
#' **Conditions.** The reference law used for this statistic treats the
#' \eqn{q} componentwise ratios as independent.  Diagonalizing the lag-zero
#' covariance does not in general diagonalize the limiting long-run
#' covariance: the independent-component law applies when \eqn{\hat L^{-1}\Omega \hat L^{-\top}} is diagonal in the limit, where
#' \eqn{\Omega} is the long-run covariance of the influence contributions.
#' The prewhitening is partial because serial and lead-lag dependence remain.
#' The statistic depends on the order of the coordinates, and it is not affine
#' equivariant: a rotation or shear of the parameters changes its value.  Use
#' `method = "hull"` when equivariance under general nonsingular
#' reparameterization is required.
#'
#' @param object An [aersn()] fit.
#' @param order Optional permutation of `1:q` giving the coordinate order used
#'   for the factorization.  The default keeps the order of the fit.  The
#'   statistic is reported in the original coordinates whichever order is used,
#'   but its value generally changes with the order.
#'
#' @return An [aersn_normalizer()] object of class
#'   `"aersn_normalizer_ldl"`.  In addition to the common fields it carries
#'   `lag_zero` (the sample lag-zero covariance in factorization order),
#'   `ldl` (the factorization in that same order),
#'   `transform` (the matrix applying the coordinate permutation followed by
#'   \eqn{\hat L^{-1}} to the error and path in their original column order), `transformed_path`, `component_ranges` and, in `diagnostics`, the
#'   reconstruction error, the maximum absolute off-diagonal element of the
#'   transformed lag-zero covariance, and the Cholesky agreement check.
#' @references
#' Hong, Y., Linton, O., McCabe, B., Sun, J. and Wang, S. (2024).
#' Kolmogorov-Smirnov type testing for structural breaks: a new
#' adjusted-range based self-normalization approach.
#' *Journal of Econometrics*, 238(2), 105603.
#' \doi{10.1016/j.jeconom.2023.105603}
#' @seealso [aersn_normalizer()], [aersn_ldl()], [aersn_compare()].
#' @examples
#' set.seed(2)
#' fit <- aersn_mean(matrix(rnorm(600), 200, 3))
#' nz <- aersn_ldl_normalizer(fit)
#' nz$ldl$L
#' nz$component_ranges
#' aersn_test(fit, method = "ldl", draws = 400, seed = 1)
#' @export
aersn_ldl_normalizer <- function(object, order = NULL) {
  check_aersn(object)
  q <- object$q
  perm <- if (is.null(order)) seq_len(q) else check_permutation(order, q)
  psi <- object$psi[, perm, drop = FALSE]
  G <- object$path$G[, perm, drop = FALSE]
  n <- object$n
  centered <- sweep(psi, 2L, colMeans(psi), "-")
  Gamma0 <- crossprod(centered) / n
  fac <- aersn_ldl(Gamma0, "the sample lag-zero covariance of the influence contributions")
  Linv <- forwardsolve(fac$L, diag(q))
  Gt <- G %*% t(Linv)
  ranges <- apply(Gt, 2L, function(x) max(x) - min(x))
  if (any(ranges <= 0)) {
    aersn_abort("Transformed coordinate(s) ",
                paste(which(ranges <= 0), collapse = ", "),
                " have zero adjusted range after LDL prewhitening; the ",
                "componentwise statistic is undefined.")
  }
  ## V in the permuted coordinates, then mapped back to the fit's order.
  Vperm <- fac$L %*% diag(ranges^2, nrow = q) %*% t(fac$L)
  inv_perm <- order(perm)
  V <- Vperm[inv_perm, inv_perm, drop = FALSE]
  dimnames(V) <- list(object$names, object$names)
  ## Diagnostics: the transformed lag-zero covariance should be diagonal;
  ## the transformed long-run covariance need not be.
  Gamma0_t <- Linv %*% Gamma0 %*% t(Linv)
  offdiag <- max(abs(Gamma0_t - diag(diag(Gamma0_t), nrow = q)))
  ## Cholesky route, for the documented scale-cancellation check.
  Cinv <- forwardsolve(fac$chol, diag(q))
  ranges_chol <- apply(G %*% t(Cinv), 2L, function(x) max(x) - min(x))
  chol_gap <- max(abs(ranges_chol * sqrt(fac$D) - ranges))
  # Columns map the ORIGINAL coordinate order into the permuted LDL path.
  transform <- Linv[, inv_perm, drop = FALSE]
  dimnames(transform) <- list(NULL, object$names)
  nz <- new_normalizer(
    "ldl", q = q, n = n, names = object$names, matrix = V,
    tuning = list(order = perm,
                  lag_zero_divisor = "n",
                  factor = "unit lower triangular L"),
    diagnostics = list(
      ldl_reconstruction_error = fac$reconstruction_error,
      transformed_lag_zero_max_offdiagonal = offdiag,
      cholesky_scale_cancellation_gap = chol_gap,
      component_range_spread = max(ranges) / min(ranges)),
    centering = object$centering, nodes = object$nodes, source = object,
    notes = c(paste0("Diagonalizing the sample lag-zero covariance does not ",
                     "imply a diagonal long-run covariance; the independent-",
                     "component reference law needs that additional condition."),
              "Not affine equivariant: the value depends on coordinate order and on nonsingular reparameterization.")
  )
  add_ldl_parts(nz, Gamma0 = Gamma0, fac = fac, transform = transform,
                transformed_path = Gt, ranges = ranges, perm = perm)
}

add_ldl_parts <- function(nz, Gamma0, fac, transform, transformed_path,
                          ranges, perm) {
  nz$lag_zero <- Gamma0
  nz$ldl <- fac
  nz$transform <- transform
  nz$transformed_path <- transformed_path
  nz$component_ranges <- stats::setNames(ranges, nz$names[perm])
  nz
}

check_permutation <- function(order, q) {
  if (!is.numeric(order) || length(order) != q || any(!is.finite(order)) ||
      any(order != floor(order))) {
    aersn_abort("`order` must be a permutation of 1:", q, ".")
  }
  o <- as.integer(order)
  if (!setequal(o, seq_len(q))) {
    aersn_abort("`order` must be a permutation of 1:", q,
                "; got ", paste(order, collapse = ", "), ".")
  }
  o
}
