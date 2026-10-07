test_that("finite-sample coverage of simultaneous intervals is near nominal (not on CRAN)", {
  skip_on_cran()
  set.seed(25)
  n <- 100; q <- 2; reps <- 300
  ref <- aersn_reference(q, n = n, draws = 3000, seed = 7)
  cv <- aersn_critical_value(ref, 0.95)
  cover_region <- logical(reps)
  cover_proj <- logical(reps)
  for (r in seq_len(reps)) {
    Y <- matrix(rnorm(n * q), n, q)
    fit <- aersn_mean(Y)
    cover_region[r] <- aersn_gauge(fit, c(0, 0)) <= cv
    ci <- confint(fit, reference = ref)
    cover_proj[r] <- all(ci[, "lower"] <= 0 & ci[, "upper"] >= 0)
  }
  ## binomial standard error about 0.013 at 300 replications
  expect_gt(mean(cover_region), 0.90)
  expect_lt(mean(cover_region), 0.99)
  expect_gte(mean(cover_proj), mean(cover_region))
})
