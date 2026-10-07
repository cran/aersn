test_that("the quadratic normalizer matches the replication implementation", {
  for (q in c(1L, 2L, 3L, 5L, 10L)) {
    Y <- sim_var1(200, q, phi = 0.4, seed = 40 + q)
    fit <- aersn_mean(Y)
    nz <- aersn_shao_normalizer(fit)
    expect_equal(unname(nz$matrix), unname(ref_shao_matrix(fit$psi)),
                 tolerance = 1e-12)
    z <- sqrt(fit$n) * fit$estimate
    expect_equal(normalizer_statistic(nz, matrix(z, nrow = 1L)),
                 ref_safe_qform(z, ref_shao_matrix(fit$psi)),
                 tolerance = 1e-10)
    expect_identical(nz$tuning$integration, "calendar")
  }
})

test_that("for one parameter the quadratic statistic is the squared ratio", {
  y <- as.numeric(sim_var1(200, 1, phi = 0.5, seed = 41))
  fit <- aersn_mean(y)
  G <- fit$path$G[-1L, , drop = FALSE]
  z <- sqrt(fit$n) * fit$estimate
  expect_equal(normalizer_statistic(aersn_shao_normalizer(fit),
                                    matrix(z, nrow = 1L)),
               unname((abs(z) / sqrt(mean(G^2)))^2), tolerance = 1e-10)
})

test_that("variance-time integration uses the profile increments", {
  set.seed(42)
  n <- 150
  s <- matrix(rnorm(n * 2), n, 2)
  s <- sweep(s, 2L, colMeans(s), "-")            # fitted scores sum to zero
  tau <- aersn_opg_profile(s)
  fit <- aersn(rep(0, 2), s, profile = tau)
  cal <- aersn_shao_normalizer(fit, integration = "calendar")
  pro <- aersn_shao_normalizer(fit, integration = "profile")
  G <- fit$path$G[-1L, , drop = FALSE]
  expect_equal(unname(cal$matrix), unname(crossprod(G) / n), tolerance = 1e-12)
  w <- diff(fit$nodes)
  expect_equal(unname(pro$matrix), unname(crossprod(G * sqrt(w))),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(cal$matrix, pro$matrix)))
  ## The default follows the fit, and a linear profile makes the two agree.
  expect_identical(aersn_shao_normalizer(fit)$tuning$integration, "profile")
  flat <- aersn(rep(0, 2), s, profile = (0:n) / n)
  expect_identical(flat$centering, "calendar")
  expect_error(aersn_shao_normalizer(flat, integration = "profile"),
               "variance-accumulation profile")
})

test_that("the reference law is keyed by the integration rule", {
  set.seed(43)
  n <- 60
  s <- sweep(matrix(rnorm(n * 2), n, 2), 2L, 0, "-")
  s <- sweep(s, 2L, colMeans(s), "-")
  tau <- as.numeric(aersn_opg_profile(s))
  fit <- aersn(rep(0, 2), s, profile = tau)
  r_pro <- aersn_reference(2, nodes = tau, draws = 100, seed = 1,
                           statistic = "shao_sq2",
                           args = list(integration = "profile"), cache = FALSE)
  r_cal <- aersn_reference(2, nodes = tau, draws = 100, seed = 1,
                           statistic = "shao_sq2",
                           args = list(integration = "calendar"), cache = FALSE)
  expect_false(isTRUE(all.equal(r_pro$draws, r_cal$draws)))
  nz <- aersn_shao_normalizer(fit, integration = "profile")
  expect_error(aersn_test(fit, method = "shao", integration = "profile",
                          reference = r_cal), "calendar-time integration")
  expect_s3_class(aersn_test(fit, method = "shao", integration = "profile",
                             reference = r_pro), "htest")
})

test_that("the fixed-b matrix matches the replication implementation", {
  for (q in c(1L, 2L, 3L, 5L)) {
    Y <- sim_var1(200, q, phi = 0.4, seed = 50 + q)
    fit <- aersn_mean(Y)
    for (b in c(0.2, 0.5, 0.7, 0.9, 1)) {
      nz <- aersn_fixed_b_normalizer(fit, b = b)
      expect_equal(unname(nz$matrix),
                   unname(ref_fixed_b_matrix(fit$path$G, b)),
                   tolerance = 1e-12)
      expect_identical(nz$tuning$m, as.integer(round(b * fit$n)))
      expect_equal(nz$tuning$b_grid, nz$tuning$m / fit$n)
      ev <- eigen(nz$matrix, symmetric = TRUE, only.values = TRUE)$values
      expect_gt(min(ev), 0)                      # Bartlett weights are psd
    }
  }
})

test_that("b is validated and the realized b is reported", {
  fit <- aersn_mean(sim_var1(100, 2, seed = 55))
  expect_error(aersn_fixed_b_normalizer(fit, b = 0), "strictly positive")
  expect_error(aersn_fixed_b_normalizer(fit, b = 1.2), "must lie in")
  expect_error(aersn_fixed_b_normalizer(fit, b = 0.001),
               "lag truncation round\\(b \\* n\\) = 0")
  ## Rounding is reported rather than hidden.
  nz <- aersn_fixed_b_normalizer(fit, b = 0.333)
  expect_identical(nz$tuning$m, 33L)
  expect_equal(nz$tuning$b_grid, 0.33)
  expect_equal(nz$diagnostics$b_rounding, 0.33 - 0.333)
  expect_output(print(nz), "b_grid")
})

test_that("fixed-b at b = 1 is exactly twice the quadratic normalizer", {
  for (q in c(1L, 2L, 3L, 5L)) {
    Y <- sim_var1(150, q, phi = 0.5, seed = 60 + q)
    fit <- aersn_mean(Y)
    Vq <- aersn_shao_normalizer(fit)$matrix
    V1 <- aersn_fixed_b_normalizer(fit, b = 1)$matrix
    expect_equal(unname(V1), unname(2 * Vq), tolerance = 1e-12)
    z <- sqrt(fit$n) * fit$estimate
    s_shao <- normalizer_statistic(aersn_shao_normalizer(fit),
                                   matrix(z, nrow = 1L))
    s_fb <- normalizer_statistic(aersn_fixed_b_normalizer(fit, b = 1),
                                 matrix(z, nrow = 1L))
    expect_equal(s_shao, 2 * s_fb, tolerance = 1e-10)
  }
})

test_that("fixed-b at b = 1 and quadratic self-normalization give the same test", {
  ## Paired reference draws: the two statistics are a fixed multiple of each
  ## other pathwise, so the decision and the p-value must agree exactly.
  ## Independent Monte Carlo runs could hide a scaling error.
  n <- 80
  q <- 2
  nodes <- (0:n) / n
  draws <- 300
  r_shao <- aersn_reference(q, n = n, draws = draws, seed = 3,
                            statistic = "shao_sq2",
                            args = list(integration = "calendar"),
                            cache = FALSE)
  r_fb <- aersn_reference(q, n = n, draws = draws, seed = 3,
                          statistic = "fixedb_bartlett",
                          args = list(b_grid = 1, m = n), cache = FALSE)
  expect_equal(r_shao$draws, 2 * r_fb$draws, tolerance = 1e-10)
  expect_equal(aersn_critical_value(r_shao, 0.95),
               2 * aersn_critical_value(r_fb, 0.95), tolerance = 1e-10)
  Y <- sim_var1(n, q, phi = 0.5, seed = 61)
  fit <- aersn_mean(Y)
  t_shao <- aersn_test(fit, method = "shao", reference = r_shao)
  t_fb <- aersn_test(fit, method = "fixedb", b = 1, reference = r_fb)
  expect_equal(unname(t_shao$statistic), unname(2 * t_fb$statistic),
               tolerance = 1e-10)
  expect_equal(t_shao$p.value, t_fb$p.value)
  expect_identical(t_shao$reject, t_fb$reject)
  ## Interval widths agree as well, since the critical values scale together.
  expect_equal(unclass(confint(fit, method = "shao", reference = r_shao))[, ],
               unclass(confint(fit, method = "fixedb", b = 1,
                               reference = r_fb))[, ], tolerance = 1e-9)
})

test_that("fixed-b and the Bartlett HAC kernel agree under their own conventions", {
  ## The two are computed by different routes: fixed-b from the partial-sum
  ## identity, HAC from the weighted sample autocovariances. They agree when
  ## the bandwidth conventions are lined up. Fixed-b uses w = 1 - lag/(b n),
  ## so a realized fraction m/n is the real bandwidth h = m, which for the
  ## Bartlett kernel is the lag truncation m - 1, not m.
  n <- 200L
  fit <- aersn_mean(sim_var1(n, 2, phi = 0.5, seed = 62))
  for (m in c(20L, 50L, 100L)) {
    Vfb <- aersn_fixed_b_normalizer(fit, b = m / n)$matrix
    expect_equal(unname(Vfb),
                 unname(aersn_hac_lrv(fit, kernel = "Bartlett",
                                      bandwidth = m)$matrix),
                 tolerance = 1e-10, info = paste("m =", m))
    expect_equal(unname(Vfb),
                 unname(aersn_hac_lrv(fit, kernel = "Bartlett",
                                      lag = m - 1L)$matrix),
                 tolerance = 1e-10)
    ## Lag truncation m is a different estimator: the conventions differ by one.
    expect_false(isTRUE(all.equal(
      unname(Vfb),
      unname(aersn_hac_lrv(fit, kernel = "Bartlett", lag = m)$matrix))))
  }
  ## The two methods still carry different reference laws.
  expect_identical(aersn_fixed_b_normalizer(fit, b = 0.25)$reference_family,
                   "fixedb_bartlett")
  expect_identical(aersn_hac_lrv(fit, bandwidth = 50)$reference_family, "chisq")
})
