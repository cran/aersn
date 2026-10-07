test_that("every method builds a normalizer with the documented contract", {
  fit <- aersn_mean(sim_var1(200, 2, seed = 1))
  expected <- list(
    hull = c(family = "geometric", scale = "gauge", ref = "hull_gauge"),
    ldl = c(family = "quadratic", scale = "wald", ref = "componentwise_mq2"),
    shao = c(family = "quadratic", scale = "wald", ref = "shao_sq2"),
    hac = c(family = "quadratic", scale = "wald", ref = "chisq"),
    fixedb = c(family = "quadratic", scale = "wald", ref = "fixedb_bartlett"),
    ewc = c(family = "quadratic", scale = "wald", ref = "ewc_F"))
  expect_identical(aersn_methods(),
                   c("hull", "ldl", "shao", "hac", "fixedb", "ewc"))
  for (m in aersn_methods()) {
    nz <- aersn_normalizer(fit, m)
    expect_s3_class(nz, "aersn_normalizer")
    expect_s3_class(nz, paste0("aersn_normalizer_", m))
    expect_identical(nz$family, unname(expected[[m]]["family"]))
    expect_identical(nz$scale, unname(expected[[m]]["scale"]))
    expect_identical(nz$reference_family, unname(expected[[m]]["ref"]))
    expect_identical(nz$q, 2L)
    expect_identical(nz$n, 200L)
    expect_identical(nz$names, fit$names)
    if (nz$family == "quadratic") {
      expect_true(is.matrix(nz$matrix))
      expect_equal(dim(nz$matrix), c(2L, 2L))
      expect_equal(nz$matrix, t(nz$matrix))
      ev <- eigen(nz$matrix, symmetric = TRUE, only.values = TRUE)$values
      expect_gt(min(ev), 0)
      expect_null(nz$hull)
    } else {
      expect_s3_class(nz$hull, "aersn_hull")
      expect_null(nz$matrix)
    }
    expect_output(print(nz), nz$label, fixed = TRUE)
  }
})

test_that("the statistic and the interval radius agree with their definitions", {
  fit <- aersn_mean(sim_var1(150, 3, seed = 2))
  z <- c(0.3, -0.2, 0.5)
  A <- rbind(c(1, 0, 0), c(1, -1, 0), c(0.2, 0.4, -0.1))
  for (m in aersn_methods()) {
    nz <- aersn_normalizer(fit, m)
    stat <- normalizer_statistic(nz, matrix(z, nrow = 1L))
    if (nz$family == "quadratic") {
      expect_equal(stat, drop(t(z) %*% solve(nz$matrix, z)))
      ## The half-width solves sup{a'u/sqrt(n) : u'V^{-1}u <= cv}.
      expect_equal(normalizer_radius(nz, A, 4),
                   sqrt(4 * rowSums((A %*% nz$matrix) * A)) / sqrt(fit$n))
    } else {
      expect_equal(stat, aersn_gauge(fit$hull, matrix(z, nrow = 1L)))
      expect_equal(normalizer_radius(nz, A, 4),
                   4 / sqrt(fit$n) * aersn_support(fit$hull, A))
    }
    ## Homogeneity: degree one for the gauge, degree two for Wald statistics.
    deg <- if (nz$scale == "gauge") 1 else 2
    expect_equal(normalizer_statistic(nz, matrix(2.5 * z, nrow = 1L)),
                 2.5^deg * stat)
    ## Central symmetry.
    expect_equal(normalizer_statistic(nz, matrix(-z, nrow = 1L)), stat)
  }
})

test_that("unknown methods and unknown tuning arguments are rejected", {
  fit <- aersn_mean(sim_var1(120, 2, seed = 3))
  expect_error(aersn_normalizer(fit, "HAC"), "Unknown method")
  expect_error(aersn_normalizer(fit, c("hac", "ewc")), "single method")
  expect_error(aersn_normalizer(fit, "hac", kernell = "Bartlett"), "unused")
  expect_error(aersn_normalizer(fit, "fixedb", nu = 4), "unused")
  expect_error(aersn_normalizer(fit, "ewc", b = 0.5), "unused")
  expect_error(aersn_test(fit, method = "shao", kernel = "Bartlett"), "unused")
  expect_error(aersn_normalizer(fit, "hac", kernel = "Gaussian"), "should be one of")
})

test_that("a prebuilt normalizer is reused and checked against the fit", {
  fit <- aersn_mean(sim_var1(120, 2, seed = 4))
  other <- aersn_mean(sim_var1(120, 3, seed = 4))
  nz <- aersn_hac_lrv(fit, kernel = "Parzen", bandwidth = 6)
  t1 <- aersn_test(fit, method = nz)
  t2 <- aersn_test(fit, method = "hac", kernel = "Parzen", bandwidth = 6)
  expect_equal(t1$statistic, t2$statistic)
  expect_error(aersn_test(other, method = nz), "q = 2")
  expect_error(aersn_test(fit, method = nz, kernel = "Bartlett"),
               "prebuilt normalizer")
})

test_that("methods defined on calendar time refuse a profile-centered fit", {
  set.seed(5)
  psi <- matrix(rnorm(200), 100, 2)
  tau <- ((0:100) / 100)^2
  fit <- aersn(colMeans(psi), psi, profile = tau)
  for (m in c("hac", "fixedb", "ewc")) {
    expect_error(aersn_normalizer(fit, m), "calendar-time grid")
  }
  expect_s3_class(aersn_normalizer(fit, "hull"), "aersn_normalizer")
  expect_s3_class(aersn_normalizer(fit, "ldl"), "aersn_normalizer")
  expect_s3_class(aersn_normalizer(fit, "shao"), "aersn_normalizer")
})

test_that("a singular normalizer matrix stops with a diagnostic and no ridge", {
  set.seed(6)
  n <- 80
  psi <- matrix(rnorm(n * 2), n, 2)
  psi <- cbind(psi, psi[, 1] + psi[, 2])       # exactly collinear third column
  fit <- suppressWarnings(aersn(rep(0, 3), psi))
  for (m in c("shao", "hac", "fixedb", "ewc")) {
    nz <- aersn_normalizer(fit, m)
    expect_error(normalizer_statistic(nz, matrix(c(1, 1, 1), nrow = 1L)),
                 "numerically singular")
    expect_error(normalizer_statistic(nz, matrix(c(1, 1, 1), nrow = 1L)),
                 "No ridge term")
  }
  expect_error(aersn_ldl_normalizer(fit), "not numerically positive definite")
})
