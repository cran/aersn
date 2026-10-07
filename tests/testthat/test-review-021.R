test_that("normalizers cannot be reused on different data or profiles", {
  set.seed(2101)
  y <- matrix(rnorm(180), 90, 2)
  f <- aersn_mean(y)
  other <- aersn_mean(y * 2)
  prof <- aersn_mean(y, profile = (0:90 / 90)^2)
  for (m in aersn_methods()) {
    nz <- aersn_normalizer(f, m)
    expect_error(aersn_test(other, method = nz, draws = 20), "different influence")
    expect_error(aersn_test(prof, method = nz, draws = 20), "different influence")
    expect_equal(aersn_test(f, method = nz, draws = 20)$statistic,
                 aersn_test(f, method = m, draws = 20)$statistic)
    expect_error(aersn_normalizer(f, nz, b = .2), "prebuilt")
  }
})

test_that("direct constructors enforce the same profile restrictions", {
  set.seed(2102)
  f <- aersn_mean(matrix(rnorm(180), 90, 2), profile = (0:90 / 90)^5)
  for (fun in list(aersn_hac_lrv, aersn_fixed_b_normalizer, aersn_ewc_lrv)) {
    expect_error(fun(f), "calendar-time")
  }
})

test_that("singular quadratic matrices cannot produce joint regions or intervals", {
  set.seed(2103)
  y <- matrix(rnorm(180), 90, 2)
  f <- suppressWarnings(aersn_mean(cbind(y, y[, 1] + y[, 2])))
  for (m in c("shao", "hac", "fixedb", "ewc")) {
    expect_error(aersn_region(f, method = m, draws = 20), "singular")
    expect_error(confint(f, method = m, draws = 20), "singular")
    expect_true(all(is.finite(confint(f, method = m, type = "marginal", draws = 20))))
  }
})

test_that("coordinate units do not create false singularity", {
  set.seed(2104)
  y <- matrix(rnorm(180), 90, 2)
  f <- aersn_mean(y)
  units <- c(1e-5, 1e5)
  g <- aersn_mean(sweep(y, 2, units, "*"))
  for (m in aersn_methods()) {
    z1 <- matrix(sqrt(f$n) * f$estimate, 1)
    z2 <- matrix(sqrt(g$n) * g$estimate, 1)
    expect_equal(normalizer_statistic(aersn_normalizer(f, m), z1),
                 normalizer_statistic(aersn_normalizer(g, m), z2), tolerance = 1e-7)
  }
})

test_that("LDL permutation is integer and its returned transform uses original columns", {
  set.seed(2105)
  f <- aersn_mean(matrix(rnorm(270), 90, 3))
  expect_error(aersn_ldl_normalizer(f, order = c(1.2, 2.8, 3)), "permutation")
  nz <- aersn_ldl_normalizer(f, order = c(3, 1, 2))
  expect_equal(unname(f$path$G %*% t(nz$transform)), unname(nz$transformed_path))
  z <- sqrt(f$n) * f$estimate
  expect_equal(unname(sum((drop(nz$transform %*% z) / nz$component_ranges)^2)),
               unname(normalizer_statistic(nz, matrix(z, 1))))
  expect_error(aersn_ldl(matrix(numeric(), 0, 0)), "square")
})

test_that("fixed-b reference metadata matches the simulated bandwidth", {
  expect_error(aersn_reference(2, n = 40, draws = 20, statistic = "fixedb_bartlett",
                               args = list(b_grid = .5, m = 8)), "b_grid = m/n")
  expect_error(aersn_reference(2, n = 40, draws = 20, statistic = "fixedb_bartlett",
                               args = list(b_grid = .5, m = 20),
                               nodes = (0:40 / 40)^2), "uniform")
  expect_error(aersn_reference(2, n = 40, statistic = "shao_sq2",
                               args = list(integration = c("calendar", "profile"))),
               "integration")
  expect_error(aersn_reference(2, n = 40, args = list(1)), "unique")
})

test_that("quadratic spectral weights are stable near zero", {
  x <- c(0, 1e-12, 1e-9, 1e-6)
  y <- 6 * pi * x / 5
  expect_equal(aersn_kernel_weight(x, "Quadratic Spectral"),
               1 - y^2 / 10, tolerance = 1e-14)
})

test_that("prewhitened HAC is stable under changes of coordinate units", {
  set.seed(2106)
  y <- matrix(rnorm(180), 90, 2)
  f <- aersn_mean(y)
  g <- aersn_mean(sweep(y, 2, c(1e-5, 1e5), "*"))
  for (kernel in c("Bartlett", "Parzen", "Quadratic Spectral")) {
    t1 <- aersn_test(f, method = "hac", prewhite = TRUE, kernel = kernel,
                    bandwidth = 5)
    t2 <- aersn_test(g, method = "hac", prewhite = TRUE, kernel = kernel,
                    bandwidth = 5)
    expect_equal(t1$statistic, t2$statistic, tolerance = 1e-7)
    expect_equal(t1$p.value, t2$p.value, tolerance = 1e-7)
  }
})

test_that("comparison settings must identify each method unambiguously", {
  set.seed(2107)
  f <- aersn_mean(rnorm(90))
  expect_error(aersn_compare(f, settings = list(hac = list(), hac = list())),
               "named list")
  expect_error(aersn_compare(f, settings = setNames(list(list()), NA_character_)),
               "named list")
})


test_that("marginal intervals record separately selected bandwidths", {
  set.seed(2108)
  y <- cbind(rnorm(120), as.numeric(arima.sim(list(ar = .8), 120)))
  f <- aersn_mean(y)
  ci <- confint(f, type = "marginal", method = "hac", bandwidth = "andrews")
  tuning <- attr(ci, "tuning_by_contrast")
  expect_length(tuning, 2L)
  for (i in 1:2) {
    expected <- aersn_hac_lrv(aersn_mean(y[, i]), bandwidth = "andrews")$tuning
    expect_equal(tuning[[i]], expected)
  }
  expect_false(isTRUE(all.equal(tuning[[1]]$bandwidth, tuning[[2]]$bandwidth)))
})


test_that("interval output identifies the method actually used", {
  set.seed(2109)
  f <- aersn_mean(matrix(rnorm(180), 90, 2))
  for (method in c("hac", "ewc", "fixedb", "shao", "ldl", "hull")) {
    ci <- confint(f, method = method, draws = 20)
    ct <- aersn_contrast(f, c(1, -1), method = method, draws = 20)
    expect_output(print(ci), method_spec(method)$label, fixed = TRUE)
    expect_output(print(ct), method_spec(method)$label, fixed = TRUE)
  }
})
