test_that("reference caches are separated by method and by tuning value", {
  aersn_clear_cache()
  keys <- c(
    aersn_reference(2, n = 40, draws = 60, seed = 1)$key,
    aersn_reference(2, n = 40, draws = 60, seed = 1,
                    statistic = "componentwise_mq2")$key,
    aersn_reference(2, n = 40, draws = 60, seed = 1, statistic = "shao_sq2",
                    args = list(integration = "calendar"))$key,
    aersn_reference(2, n = 40, draws = 60, seed = 1, statistic = "shao_sq2",
                    args = list(integration = "profile"))$key,
    aersn_reference(2, n = 40, draws = 60, seed = 1,
                    statistic = "fixedb_bartlett",
                    args = list(b_grid = 0.5, m = 20))$key,
    aersn_reference(2, n = 40, draws = 60, seed = 1,
                    statistic = "fixedb_bartlett",
                    args = list(b_grid = 0.25, m = 10))$key)
  expect_identical(length(unique(keys)), 6L)
  expect_gte(aersn_clear_cache(), 6L)
  ## A repeated call returns the cached object unchanged.
  r1 <- aersn_reference(2, n = 40, draws = 60, seed = 1, statistic = "shao_sq2",
                        args = list(integration = "calendar"))
  r2 <- aersn_reference(2, n = 40, draws = 60, seed = 1, statistic = "shao_sq2",
                        args = list(integration = "calendar"))
  expect_identical(r1, r2)
  aersn_clear_cache()
})

test_that("reference arguments are validated", {
  expect_error(aersn_reference(2, n = 30, draws = 10, statistic = "wald"),
               "`statistic` must be one of")
  expect_error(aersn_reference(2, n = 30, draws = 10, statistic = "shao_sq2"),
               "needs reference argument")
  expect_error(aersn_reference(2, n = 30, draws = 10,
                               statistic = "componentwise_mq2",
                               args = list(b_grid = 1)), "do not apply")
  expect_error(aersn_reference(2, n = 30, draws = 10, statistic = "shao_sq2",
                               args = list(integration = "clock")),
               "calendar")
  expect_error(aersn_reference(2, n = 30, draws = 10,
                               statistic = "fixedb_bartlett",
                               args = list(b_grid = 0.5, m = 0)), "`m`")
  expect_error(aersn_reference(2, n = 30, draws = 10,
                               statistic = "shao_sq2", grid = "continuous",
                               args = list(integration = "calendar")),
               "increment-hull gauge")
})

test_that("a reference for one method is refused by another", {
  fit <- aersn_mean(sim_var1(60, 2, seed = 120))
  r_hull <- aersn_reference(2, n = 60, draws = 80, seed = 1)
  r_shao <- aersn_reference(2, n = 60, draws = 80, seed = 1,
                            statistic = "shao_sq2",
                            args = list(integration = "calendar"))
  r_fb5 <- aersn_reference(2, n = 60, draws = 80, seed = 1,
                           statistic = "fixedb_bartlett",
                           args = list(b_grid = 0.5, m = 30))
  expect_error(aersn_test(fit, method = "shao", reference = r_hull),
               "not interchangeable")
  expect_error(aersn_test(fit, method = "hull", reference = r_shao),
               "not interchangeable")
  expect_error(aersn_test(fit, method = "fixedb", b = 0.2, reference = r_fb5),
               "b = 0.5")
  expect_s3_class(aersn_test(fit, method = "fixedb", b = 0.5,
                             reference = r_fb5), "htest")
  ## Closed-form and tabulated references apply to the hull only.
  expect_error(aersn_test(fit, method = "shao", reference = "continuous"),
               "increment-hull gauge")
  expect_error(aersn_test(fit, method = "ldl", reference = "registry"),
               "registry tabulates")
  ## HAC and EWC use closed-form laws and refuse a simulated one.
  expect_error(aersn_test(fit, method = "hac", reference = r_hull),
               "does not match")
  expect_error(aersn_test(fit, method = "ewc", reference = "registry"),
               "does not apply")
})

test_that("the EWC reference is tied to its number of cosine terms", {
  fit <- aersn_mean(sim_var1(200, 2, seed = 121))
  ref16 <- aersn_parametric_reference("ewc_F", df = 2, nu = 16)
  expect_s3_class(aersn_test(fit, method = "ewc", nu = 16, reference = ref16),
                  "htest")
  expect_error(aersn_test(fit, method = "ewc", nu = 20, reference = ref16),
               "does not match")
  ## The HAC reference is tied to the dimension.
  ref3 <- aersn_parametric_reference("chisq", df = 3)
  expect_error(aersn_test(fit, method = "hac", reference = ref3),
               "does not match")
})

test_that("simulated laws are reproducible and leave the user's stream alone", {
  for (spec in list(list(s = "componentwise_mq2", a = list()),
                    list(s = "shao_sq2", a = list(integration = "calendar")),
                    list(s = "fixedb_bartlett",
                         a = list(b_grid = 0.5, m = 25)))) {
    r1 <- aersn_reference(2, n = 50, draws = 80, seed = 9,
                          statistic = spec$s, args = spec$a, cache = FALSE)
    r2 <- aersn_reference(2, n = 50, draws = 80, seed = 9,
                          statistic = spec$s, args = spec$a, cache = FALSE)
    expect_identical(r1$draws, r2$draws, info = spec$s)
    r3 <- aersn_reference(2, n = 50, draws = 80, seed = 10,
                          statistic = spec$s, args = spec$a, cache = FALSE)
    expect_false(identical(r1$draws, r3$draws))
    if (.Platform$OS.type == "unix") {
      r4 <- aersn_reference(2, n = 50, draws = 80, seed = 9,
                            statistic = spec$s, args = spec$a, cache = FALSE,
                            cores = 2)
      expect_identical(r1$draws, r4$draws)
    }
  }
  set.seed(4242)
  before <- stats::rnorm(3)
  set.seed(4242)
  invisible(aersn_reference(2, n = 40, draws = 60, seed = 1,
                            statistic = "shao_sq2",
                            args = list(integration = "calendar"),
                            cache = FALSE))
  expect_identical(before, stats::rnorm(3))
})

test_that("a simulated law with failed draws is not returned", {
  ## A grid with too few positive increments makes the bridge degenerate.
  expect_error(aersn_reference(3, n = 3, draws = 5,
                               statistic = "componentwise_mq2"),
               "needs n >= q \\+ 1 intervals")
  ## Enough grid intervals, but only two of them carry positive variance.
  nodes <- c(0, rep(0.5, 3), 1)
  expect_error(aersn_reference(2, nodes = nodes, draws = 5,
                               statistic = "shao_sq2",
                               args = list(integration = "calendar")),
               "q \\+ 1 positive")
})

test_that("printing a reference names its statistic", {
  r <- aersn_reference(2, n = 40, draws = 60, seed = 1, statistic = "shao_sq2",
                       args = list(integration = "calendar"), cache = FALSE)
  expect_output(print(r), "quadratic self-normalization")
  expect_output(print(r), "calendar-time integration")
  rf <- aersn_reference(2, n = 40, draws = 60, seed = 1,
                        statistic = "fixedb_bartlett",
                        args = list(b_grid = 0.5, m = 20), cache = FALSE)
  expect_output(print(rf), "Bartlett fixed-b")
  expect_output(print(rf), "b = 0.5")
})

test_that("the number of cosine terms is not captured by partial matching", {
  fit <- aersn_mean(sim_var1(200, 2, seed = 130))
  ## `nu` is a prefix of `null`: without an exact formal, partial matching
  ## would silently set the null value to the number of cosine terms.
  tt <- aersn_test(fit, method = "ewc", nu = 16)
  expect_identical(tt$tuning$nu, 16L)
  expect_equal(unname(tt$null.value), c(0, 0))
  tt2 <- aersn_test(fit, c(0.1, -0.2), method = "ewc", nu = 9)
  expect_identical(tt2$tuning$nu, 9L)
  expect_equal(unname(tt2$null.value), c(0.1, -0.2))
  expect_identical(attr(confint(fit, method = "ewc", nu = 11),
                        "tuning")$nu, 11L)
  expect_identical(aersn_region(fit, method = "ewc", nu = 12,
                                )$normalizer$tuning$nu, 12L)
  expect_identical(attr(aersn_contrast(fit, c(1, -1), method = "ewc", nu = 13),
                        "tuning")$nu, 13L)
})
