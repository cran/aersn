test_that("the kernel weights match sandwich::kweights", {
  x <- c(0, 0.1, 0.25, 0.5, 0.75, 0.99, 1, 1.5, 3)
  for (k in c("Bartlett", "Parzen", "Quadratic Spectral")) {
    expect_equal(aersn_kernel_weight(x, k), sandwich::kweights(x, kernel = k),
                 tolerance = 1e-12, info = k)
    expect_equal(aersn_kernel_weight(-x, k), aersn_kernel_weight(x, k))
  }
  ## Compact support for Bartlett and Parzen, unbounded support for QS.
  expect_equal(aersn_kernel_weight(c(1, 2, 5), "Bartlett"), rep(0, 3))
  expect_equal(aersn_kernel_weight(c(1, 2, 5), "Parzen"), rep(0, 3))
  expect_true(any(aersn_kernel_weight(c(2, 5, 12), "Quadratic Spectral") != 0))
  expect_equal(aersn_kernel_weight(0, "Quadratic Spectral"), 1)
})

test_that("the HAC matrix matches a direct implementation of the definition", {
  for (q in c(1L, 2L, 3L)) {
    Y <- sim_var1(60, q, phi = 0.5, seed = 70 + q)
    fit <- aersn_mean(Y)
    for (k in c("Bartlett", "Parzen", "Quadratic Spectral")) {
      h <- 4.5
      nz <- aersn_hac_lrv(fit, kernel = k, bandwidth = h)
      w <- aersn_kernel_weight(seq_len(nrow(Y) - 1L) / h, k)
      expect_equal(unname(nz$matrix), unname(ref_hac_matrix(fit$psi, w)),
                   tolerance = 1e-10, info = paste(k, q))
    }
  }
})

test_that("the HAC covariance and Wald statistic agree with sandwich::vcovHAC", {
  set.seed(71)
  n <- 200
  x1 <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  x2 <- stats::rnorm(n)
  y <- 1 + 0.7 * x1 - 0.3 * x2 +
    as.numeric(stats::arima.sim(list(ar = 0.6), n))
  lmfit <- stats::lm(y ~ x1 + x2)
  fit <- aersn_lm(y ~ x1 + x2)
  b0 <- c(1, 0.7, -0.3)
  for (k in c("Bartlett", "Parzen", "Quadratic Spectral")) {
    for (prewhite in c(FALSE, TRUE)) {
      nz <- aersn_hac_lrv(fit, kernel = k, bandwidth = 5, prewhite = prewhite)
      w <- sandwich::kweights(0:(n - 1L - as.integer(prewhite)) / 5, kernel = k)
      Vs <- sandwich::vcovHAC(lmfit, weights = w,
                              prewhite = as.integer(prewhite),
                              adjust = FALSE, sandwich = TRUE)
      ## The package normalizes the influence contributions, so its matrix is
      ## n times the sandwich covariance of the coefficients.
      expect_equal(unname(nz$matrix), unname(n * Vs), tolerance = 1e-9,
                   info = paste(k, prewhite))
      tt <- aersn_test(fit, null = b0, method = "hac", kernel = k,
                       bandwidth = 5, prewhite = prewhite)
      wald <- drop(t(coef(lmfit) - b0) %*% solve(Vs) %*% (coef(lmfit) - b0))
      expect_equal(unname(tt$statistic), wald, tolerance = 1e-8)
      expect_equal(tt$p.value, stats::pchisq(wald, 3, lower.tail = FALSE),
                   tolerance = 1e-10)
      expect_equal(tt$critical.value, stats::qchisq(0.95, 3))
    }
  }
  ## The finite-sample multiplier is exactly n/(n-q).
  v0 <- aersn_hac_lrv(fit, bandwidth = 5)$matrix
  v1 <- aersn_hac_lrv(fit, bandwidth = 5, adjust = TRUE)$matrix
  expect_equal(v1, v0 * n / (n - 3))
})

test_that("bandwidth rules return the documented values", {
  n <- 400
  Y <- sim_var1(n, 2, phi = 0.5, seed = 72)
  fit <- aersn_mean(Y)
  short <- aersn_hac_lrv(fit, bandwidth = "short")
  expect_identical(short$tuning$lag_truncation,
                   as.integer(floor(4 * (n / 100)^(2 / 9))))
  expect_equal(short$tuning$bandwidth, floor(4 * (n / 100)^(2 / 9)) + 1)
  long <- aersn_hac_lrv(fit, bandwidth = "long")
  expect_identical(long$tuning$lag_truncation, as.integer(ceiling(1.3 * sqrt(n))))
  expect_equal(long$tuning$bandwidth, ceiling(1.3 * sqrt(n)) + 1)
  E <- sweep(fit$psi, 2L, colMeans(fit$psi), "-")
  for (k in c("Bartlett", "Parzen", "Quadratic Spectral")) {
    a <- aersn_hac_lrv(fit, kernel = k, bandwidth = "andrews")
    expect_equal(a$tuning$bandwidth,
                 sandwich::bwAndrews(E, kernel = k, weights = rep(1, 2),
                                     prewhite = FALSE),
                 tolerance = 1e-10, info = k)
    expect_match(a$tuning$bandwidth_rule, "Andrews")
  }
  nw <- aersn_hac_lrv(fit, kernel = "Bartlett", bandwidth = "newey-west")
  expect_identical(nw$tuning$lag_truncation,
                   as.integer(floor(sandwich::bwNeweyWest(E, kernel = "Bartlett",
                                                          weights = rep(1, 2),
                                                          prewhite = FALSE))))
  expect_equal(nw$tuning$bandwidth, nw$tuning$lag_truncation + 1)
  ## Distinct rules really do give distinct bandwidths here.
  expect_false(short$tuning$bandwidth == long$tuning$bandwidth)
  expect_false(short$tuning$bandwidth == nw$tuning$bandwidth)
  ## An explicit lag truncation is h = L + 1.
  lg <- aersn_hac_lrv(fit, kernel = "Bartlett", lag = 7)
  expect_equal(lg$tuning$bandwidth, 8)
  expect_identical(lg$tuning$lag_truncation, 7L)
  expect_equal(lg$matrix, aersn_hac_lrv(fit, bandwidth = 8)$matrix)
})

test_that("the Andrews Parzen bandwidth is the quadratic spectral one rescaled", {
  ## The replication code derives the Parzen bandwidth from the quadratic
  ## spectral one through the ratio of Andrews' kernel constants.
  for (seed in 73:75) {
    Y <- sim_var1(300, 3, phi = 0.6, seed = seed)
    E <- sweep(Y, 2L, colMeans(Y), "-")
    aP <- sandwich::bwAndrews(E, kernel = "Parzen", weights = rep(1, 3),
                              prewhite = FALSE)
    aQ <- sandwich::bwAndrews(E, kernel = "Quadratic Spectral",
                              weights = rep(1, 3), prewhite = FALSE)
    expect_equal(aP, aQ * 2.6614 / 1.3221, tolerance = 1e-10)
  }
})

test_that("unsupported kernel and rule combinations are refused explicitly", {
  fit <- aersn_mean(sim_var1(200, 2, seed = 76))
  for (k in c("Parzen", "Quadratic Spectral")) {
    expect_error(aersn_hac_lrv(fit, kernel = k, bandwidth = "newey-west"),
                 "Bartlett kernel only")
  }
  expect_error(aersn_hac_lrv(fit, kernel = "Quadratic Spectral", lag = 5),
               "not defined for the")
  expect_error(aersn_hac_lrv(fit, bandwidth = "plugin"), "should be one of")
  expect_error(aersn_hac_lrv(fit, bandwidth = -1), "must lie in")
  expect_error(aersn_hac_lrv(fit, bandwidth = 0), "must be positive")
  expect_error(aersn_hac_lrv(fit, lag = 2.5), "whole number")
  expect_error(aersn_hac_lrv(fit, prewhite = "yes"), "TRUE or FALSE")
  ## Prewhitening needs enough observations to fit the VAR(1).
  small <- aersn_mean(sim_var1(4, 3, seed = 77))
  expect_error(aersn_hac_lrv(small, prewhite = TRUE), "prewhitening")
})

test_that("centering has no effect when the contributions already sum to zero", {
  fit <- aersn_mean(sim_var1(150, 2, seed = 78))
  expect_equal(max(abs(colSums(fit$psi))), 0, tolerance = 1e-10)
  expect_equal(aersn_hac_lrv(fit, bandwidth = 5, center = TRUE)$matrix,
               aersn_hac_lrv(fit, bandwidth = 5, center = FALSE)$matrix,
               tolerance = 1e-10)
  ## With contributions that do not sum to zero it does matter.
  psi <- sim_var1(150, 2, seed = 79) + 3
  fit2 <- aersn(c(0, 0), psi)
  expect_false(isTRUE(all.equal(
    aersn_hac_lrv(fit2, bandwidth = 5, center = TRUE)$matrix,
    aersn_hac_lrv(fit2, bandwidth = 5, center = FALSE)$matrix)))
})
