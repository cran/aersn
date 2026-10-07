## Independent reference implementations of the comparator normalizers,
## transcribed from the manuscript's replication code so that the package
## implementation is compared with a separately written algorithm rather than
## with itself.  Sources:
##   replication/core_methods/R/multivariate_mean_comparison.R
##     (centered_bridge, denominator_objects, safe_qform, lower_sqrt)
##   replication/core_methods/R/multivariate_fixedb_har_comparison.R
##     (fixed_b_matrix, ewc_matrix, critical_values)
##   replication/revisions/code/table3_kernel_bandwidth.R
##     (kernel weights, bandwidth rules)

ref_centered_bridge <- function(Y) {
  n <- nrow(Y)
  q <- ncol(Y)
  cs <- apply(Y, 2L, cumsum)
  if (q == 1L) cs <- matrix(cs, ncol = 1L)
  rbind(rep(0, q), cs - tcrossprod(seq_len(n) / n, cs[n, ])) / sqrt(n)
}

ref_lower_sqrt <- function(S) t(chol((S + t(S)) / 2))

ref_safe_qform <- function(z, V) {
  V <- (V + t(V)) / 2
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(ev)) return(NA_real_)
  drop(crossprod(z, solve(V, z)))
}

## Lag-zero componentwise statistic exactly as in the replication code:
## the lower Cholesky factor of the sample lag-zero covariance is used.
ref_componentwise <- function(psi, z) {
  n <- nrow(psi)
  G <- ref_centered_bridge(psi)
  V0 <- crossprod(sweep(psi, 2L, colMeans(psi), "-")) / n
  C0 <- ref_lower_sqrt(V0)
  z0 <- as.numeric(solve(C0, z))
  G0 <- G %*% solve(t(C0))
  ranges0 <- apply(G0, 2L, function(x) diff(range(x)))
  list(statistic = sum((z0 / ranges0)^2), C = C0, ranges = ranges0,
       halfwidth_matrix = C0 %*% diag(ranges0^2, nrow = ncol(psi)) %*% t(C0))
}

ref_shao_matrix <- function(psi) {
  n <- nrow(psi)
  G <- ref_centered_bridge(psi)
  crossprod(G[-1L, , drop = FALSE]) / n
}

ref_fixed_b_matrix <- function(G, b = .5) {
  N <- nrow(G) - 1L
  m <- as.integer(round(b * N))
  b_grid <- m / N
  level <- crossprod(G[-1L, , drop = FALSE]) / N
  left <- G[seq_len(N - m + 1L), , drop = FALSE]
  right <- G[(m + 1L):(N + 1L), , drop = FALSE]
  lagged <- crossprod(left, right) / N
  (2 / b_grid) * level - (lagged + t(lagged)) / b_grid
}

ref_ewc_matrix <- function(Y, nu) {
  n <- nrow(Y)
  tt <- seq_len(n) - .5
  Phi <- outer(tt, seq_len(nu), function(t, j) sqrt(2 / n) * cos(pi * j * t / n))
  Xi <- crossprod(Phi, Y)
  crossprod(Xi) / nu
}

ref_ewc_critical <- function(q, nu, level = .95) {
  nu * q / (nu - q + 1) * stats::qf(level, q, nu - q + 1)
}

## Naive double-loop HAC, written directly from the definition.
ref_hac_matrix <- function(E, w) {
  n <- nrow(E)
  q <- ncol(E)
  E <- sweep(E, 2L, colMeans(E), "-")
  V <- matrix(0, q, q)
  for (s in seq_len(n)) {
    for (t in seq_len(n)) {
      lag <- abs(s - t)
      wt <- if (lag == 0L) 1 else if (lag <= length(w)) w[lag] else 0
      if (wt != 0) V <- V + wt * tcrossprod(E[s, ], E[t, ])
    }
  }
  V / n
}

## A bivariate autoregression used by several tests.
sim_var1 <- function(n, q = 2L, phi = 0.4, seed = 1L, burn = 200L) {
  set.seed(seed)
  e <- matrix(stats::rnorm((n + burn) * q), n + burn, q)
  Y <- e
  for (t in 2:(n + burn)) Y[t, ] <- phi * Y[t - 1L, ] + e[t, ]
  Y[(burn + 1L):(burn + n), , drop = FALSE]
}
