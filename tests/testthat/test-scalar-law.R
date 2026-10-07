test_that("series and Bessel representations of the scalar law agree", {
  cs <- c(0.03, 0.1, 0.3, 0.7, 1, 1.5, 2, 2.5, 3.5)
  a <- aersn_scalar_cdf(cs, method = "series")
  b <- aersn_scalar_cdf(cs, method = "bessel")
  expect_equal(a, b, tolerance = 1e-9)
  d1 <- aersn_scalar_density(cs, method = "series")
  d2 <- aersn_scalar_density(cs, method = "bessel")
  expect_equal(d1, d2, tolerance = 1e-8)
  expect_equal(aersn_scalar_density(-cs), aersn_scalar_density(cs))
})

test_that("the density is the derivative of the |M| law and integrates to one", {
  for (c0 in c(0.4, 1.2, 2.1)) {
    num <- (aersn_scalar_cdf(c0 + 1e-5) - aersn_scalar_cdf(c0 - 1e-5)) / 2e-5
    expect_equal(num, 2 * aersn_scalar_density(c0), tolerance = 1e-6)
  }
  expect_equal(integrate(aersn_scalar_density, -Inf, Inf)$value, 1,
               tolerance = 1e-7)
  expect_equal(aersn_scalar_density(0), 0.5, tolerance = 1e-9)
  expect_equal(aersn_scalar_cdf(0), 0)
  expect_true(all(diff(aersn_scalar_cdf(seq(0, 6, by = 0.1))) > 0))
})

test_that("closed-form quantiles reproduce the manuscript's analytic table", {
  ## Supplement Table "Continuous-path scalar critical values" and the
  ## replication file scalar_M_analytic_critical_values.csv.
  target <- c(`0.9` = 1.39739012671653, `0.95` = 1.70577645482396,
              `0.975` = 1.99784763830686, `0.99` = 2.36737759542925,
              `0.995` = 2.63809633080277, `0.999` = 3.2470201381194)
  got <- aersn_scalar_quantile(as.numeric(names(target)))
  expect_equal(unname(got), unname(target), tolerance = 1e-8)
  expect_equal(aersn_scalar_cdf(got), as.numeric(names(target)),
               tolerance = 1e-10)
})

test_that("scalar law arguments are validated", {
  expect_error(aersn_scalar_cdf(-1), "nonnegative")
  expect_error(aersn_scalar_quantile(1), "strictly between")
  expect_error(aersn_scalar_density(NA), "finite")
})
