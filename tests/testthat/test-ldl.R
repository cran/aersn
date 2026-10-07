test_that("the LDL factorization reconstructs the matrix and matches Cholesky", {
  set.seed(11)
  for (q in c(1L, 2L, 3L, 5L)) {
    M <- matrix(rnorm(q * q), q, q)
    S <- crossprod(M) + diag(q)
    f <- aersn_ldl(S)
    expect_equal(f$L %*% diag(f$D, nrow = q) %*% t(f$L), S)
    expect_lt(f$reconstruction_error, 1e-10)
    expect_equal(diag(f$L), rep(1, q))            # unit lower triangular
    expect_true(all(f$L[upper.tri(f$L)] == 0))
    expect_true(all(f$D > 0))
    expect_equal(f$chol, f$L %*% diag(sqrt(f$D), nrow = q))
    expect_equal(f$chol, t(chol(S)))
  }
  expect_error(aersn_ldl(matrix(1:6, 2, 3)), "must be square")
  expect_error(aersn_ldl(matrix(c(1, 2, 3, 4), 2, 2)), "not symmetric")
  expect_error(aersn_ldl(matrix(c(1, 1, 1, 1), 2, 2)),
               "not numerically positive definite")
  expect_error(aersn_ldl(matrix(c(1, 1, 1, 1), 2, 2)), "No ridge term")
})

test_that("the componentwise statistic matches the replication implementation", {
  for (q in c(2L, 3L, 5L)) {
    Y <- sim_var1(200, q, phi = 0.5, seed = 20 + q)
    fit <- aersn_mean(Y)
    nz <- aersn_ldl_normalizer(fit)
    z <- sqrt(fit$n) * fit$estimate
    ref <- ref_componentwise(fit$psi, z)
    expect_equal(normalizer_statistic(nz, matrix(z, nrow = 1L)),
                 ref$statistic, tolerance = 1e-10)
    ## The Wald matrix is the one the replication code uses for half-widths.
    expect_equal(unname(nz$matrix), unname(ref$halfwidth_matrix),
                 tolerance = 1e-10)
    ## Unit-lower and Cholesky factors give the same statistic: the diagonal
    ## scale cancels between numerator and denominator.
    expect_equal(nz$ldl$chol, ref$C, tolerance = 1e-12)
    expect_lt(nz$diagnostics$cholesky_scale_cancellation_gap, 1e-9)
    zc <- as.numeric(solve(ref$C, z))
    expect_equal(sum((zc / ref$ranges)^2),
                 sum((as.numeric(nz$transform %*% z) / nz$component_ranges)^2),
                 tolerance = 1e-10)
  }
})

test_that("the transformation diagonalizes the lag-zero covariance only", {
  Y <- sim_var1(300, 3, phi = 0.6, seed = 31)
  fit <- aersn_mean(Y)
  nz <- aersn_ldl_normalizer(fit)
  expect_lt(nz$diagnostics$transformed_lag_zero_max_offdiagonal, 1e-10)
  ## The transformed long-run covariance is not diagonal in general, which is
  ## the condition the independent-component reference law needs.
  lrv <- aersn_hac_lrv(fit, bandwidth = "long")$matrix
  tl <- nz$transform %*% lrv %*% t(nz$transform)
  expect_gt(max(abs(tl[upper.tri(tl)])) / max(abs(diag(tl))), 1e-3)
  expect_true(any(grepl("long-run covariance", nz$notes)))
})

test_that("the componentwise statistic depends on coordinate order", {
  Y <- sim_var1(200, 3, phi = 0.5, seed = 32)
  fit <- aersn_mean(Y)
  z <- sqrt(fit$n) * fit$estimate
  s1 <- normalizer_statistic(aersn_ldl_normalizer(fit), matrix(z, nrow = 1L))
  s2 <- normalizer_statistic(aersn_ldl_normalizer(fit, order = c(3, 2, 1)),
                             matrix(z, nrow = 1L))
  expect_false(isTRUE(all.equal(s1, s2)))
  ## Reordering is a relabelling: the statistic equals the one obtained by
  ## permuting the data and using the natural order.
  fitp <- aersn_mean(Y[, c(3, 2, 1)])
  zp <- sqrt(fitp$n) * fitp$estimate
  expect_equal(normalizer_statistic(aersn_ldl_normalizer(fitp),
                                    matrix(zp, nrow = 1L)), s2,
               tolerance = 1e-10)
  expect_error(aersn_ldl_normalizer(fit, order = c(1, 1, 2)), "permutation")
  expect_error(aersn_ldl_normalizer(fit, order = 1:2), "permutation")
})

test_that("for one parameter the componentwise statistic is the squared hull statistic", {
  y <- as.numeric(sim_var1(200, 1, phi = 0.5, seed = 33))
  fit <- aersn_mean(y)
  z <- matrix(sqrt(fit$n) * fit$estimate, nrow = 1L)
  hull <- normalizer_statistic(aersn_hull_normalizer(fit), z)
  ldl <- normalizer_statistic(aersn_ldl_normalizer(fit), z)
  expect_equal(ldl, hull^2, tolerance = 1e-12)
  expect_equal(unname(aersn_ldl_normalizer(fit)$ldl$L), matrix(1, 1, 1))
  ## The two reference laws are simulated from the same draws, so their
  ## quantiles satisfy the same relation exactly.
  rh <- aersn_reference(1, n = fit$n, draws = 500, seed = 7,
                        statistic = "hull_gauge", cache = FALSE)
  rl <- aersn_reference(1, n = fit$n, draws = 500, seed = 7,
                        statistic = "componentwise_mq2", cache = FALSE)
  expect_equal(rl$draws, rh$draws^2)
  ## The p-value is exactly the same: the reference draws are a monotone
  ## transformation of each other and so is the statistic.
  expect_equal(aersn_pvalue(rl, ldl)$p.value, aersn_pvalue(rh, hull)$p.value)
  ## The critical values agree up to the linear interpolation of the type-8
  ## sample quantile, which does not commute with squaring.
  expect_equal(aersn_critical_value(rl, 0.95),
               aersn_critical_value(rh, 0.95)^2, tolerance = 1e-4)
})

test_that("the componentwise statistic does not depend on the lag-zero divisor", {
  ## The unit lower triangular factor is unchanged when the lag-zero
  ## covariance is multiplied by a positive scalar, so dividing by n or by
  ## n - 1, as two of the research scripts do, gives the same statistic.
  fit <- aersn_mean(sim_var1(200, 3, phi = 0.4, seed = 34))
  psi <- fit$psi
  z <- sqrt(fit$n) * fit$estimate
  centered <- sweep(psi, 2L, colMeans(psi), "-")
  S_n <- crossprod(centered) / fit$n
  S_n1 <- stats::cov(psi)
  S_scaled <- 7.3 * S_n
  expect_equal(aersn_ldl(S_n)$L, aersn_ldl(S_n1)$L, tolerance = 1e-12)
  expect_equal(aersn_ldl(S_n)$L, aersn_ldl(S_scaled)$L, tolerance = 1e-12)
  expect_equal(aersn_ldl(S_scaled)$D, 7.3 * aersn_ldl(S_n)$D, tolerance = 1e-10)
  stat_from <- function(S) {
    C <- t(chol(S))
    G0 <- fit$path$G %*% solve(t(C))
    sum((as.numeric(solve(C, z)) /
           apply(G0, 2L, function(x) diff(range(x))))^2)
  }
  target <- normalizer_statistic(aersn_ldl_normalizer(fit), matrix(z, nrow = 1L))
  expect_equal(stat_from(S_n), target, tolerance = 1e-10)
  expect_equal(stat_from(S_n1), target, tolerance = 1e-10)
  expect_equal(stat_from(S_scaled), target, tolerance = 1e-10)
})
