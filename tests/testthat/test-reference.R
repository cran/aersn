test_that("reference simulation is reproducible and independent of cores", {
  r1 <- aersn_reference(2, n = 30, draws = 250, seed = 11, cache = FALSE)
  r2 <- aersn_reference(2, n = 30, draws = 250, seed = 11, cache = FALSE)
  expect_identical(r1$draws, r2$draws)
  if (.Platform$OS.type == "unix") {
    r3 <- aersn_reference(2, n = 30, draws = 250, seed = 11, cache = FALSE,
                          cores = 2)
    expect_identical(r1$draws, r3$draws)
  }
  r4 <- aersn_reference(2, n = 30, draws = 250, seed = 12, cache = FALSE)
  expect_false(identical(r1$draws, r4$draws))
  s1 <- aersn_reference(1, n = 30, draws = 3000, seed = 11, cache = FALSE)
  s2 <- aersn_reference(1, n = 30, draws = 3000, seed = 11, cache = FALSE)
  expect_identical(s1$draws, s2$draws)
  expect_length(s1$signed, 3000L)
  expect_equal(abs(s1$signed), s1$draws)
})

test_that("the simulation does not disturb the user's random-number stream", {
  set.seed(123)
  a <- rnorm(3)
  set.seed(123)
  invisible(aersn_reference(2, n = 20, draws = 100, seed = 5, cache = FALSE))
  b <- rnorm(3)
  expect_identical(a, b)
})

test_that("the cache distinguishes every setting that changes the law", {
  aersn_clear_cache()
  r1 <- aersn_reference(2, n = 30, draws = 200, seed = 1)
  r2 <- aersn_reference(2, n = 30, draws = 200, seed = 1)
  expect_identical(r1, r2)
  keys <- c(aersn_reference(2, n = 30, draws = 300, seed = 1)$key,
            aersn_reference(2, n = 30, draws = 200, seed = 2)$key,
            aersn_reference(2, n = 31, draws = 200, seed = 1)$key,
            aersn_reference(3, n = 30, draws = 200, seed = 1)$key,
            aersn_reference(2, nodes = ((0:30) / 30)^2, draws = 200,
                            seed = 1)$key)
  expect_equal(length(unique(c(r1$key, keys))), 6L)
  expect_gte(aersn_clear_cache(), 6L)
})

test_that("reference objects refuse mismatched statistics", {
  set.seed(13)
  fit2 <- aersn_mean(matrix(rnorm(80), 40, 2))
  fit3 <- aersn_mean(matrix(rnorm(120), 40, 3))
  fit2b <- aersn_mean(matrix(rnorm(100), 50, 2))
  ref <- aersn_reference(2, n = 40, draws = 200, seed = 1)
  expect_s3_class(aersn_test(fit2, reference = ref), "htest")
  expect_error(aersn_test(fit3, reference = ref), "q = 2")
  expect_error(aersn_test(fit2b, reference = ref), "n = 40")
  expect_error(aersn_test(fit2, reference = "continuous"), "only for q = 1")
  ## profile-centered statistics need the profile grid
  fitp <- aersn(colMeans(matrix(rnorm(80), 40, 2)) + 0,
                matrix(rnorm(80), 40, 2), profile = ((0:40) / 40)^2)
  expect_error(aersn_test(fitp, reference = ref), "profile")
  refp <- aersn_reference(fitp, draws = 200, seed = 1)
  expect_false(refp$uniform)
  expect_s3_class(aersn_test(fitp, reference = refp), "htest")
  expect_error(aersn_test(fit2, reference = refp), "nonuniform")
  expect_error(aersn_reference(2, n = 2), "n >= q \\+ 1")
})

test_that("matched-grid quantiles reproduce the verified registry", {
  ## q = 1: registry uses 2,000,000 draws (batch MCSE about 0.0013).
  reg <- subset(aersn_registry, q == 1 & n == 200 & prob == 0.95)
  ref <- aersn_reference(1, n = 200, draws = 20000, seed = 20260911,
                         cache = FALSE)
  got <- unname(quantile(ref, 0.95))
  se <- sqrt(reg$mcse^2 + unname(aersn_mcse(ref, 0.95))^2)
  expect_lt(abs(got - reg$cv), 4 * se + 0.005)
  ## continuous-path quantile is smaller than the matched-grid quantile
  expect_lt(aersn_scalar_quantile(0.95), got)
  ## q = 2, n = 200: registry 2.5664 (MCSE 0.0129) from 30,000 draws.
  reg2 <- subset(aersn_registry, q == 2 & n == 200 & prob == 0.95)
  ref2 <- aersn_reference(2, n = 200, draws = 2000, seed = 20260911,
                          cache = FALSE)
  got2 <- unname(quantile(ref2, 0.95))
  se2 <- sqrt(reg2$mcse^2 + unname(aersn_mcse(ref2, 0.95))^2)
  expect_lt(abs(got2 - reg2$cv), 4 * se2 + 0.02)
})

test_that("larger registry cells are reproduced (not run on CRAN)", {
  skip_on_cran()
  reg <- subset(aersn_registry, q == 3 & n == 500 & prob == 0.95)
  ref <- aersn_reference(3, n = 500, draws = 3000, seed = 20260912,
                         cache = FALSE)
  got <- unname(quantile(ref, 0.95))
  se <- sqrt(reg$mcse^2 + unname(aersn_mcse(ref, 0.95))^2)
  expect_lt(abs(got - reg$cv), 4 * se + 0.02)
  reg1 <- subset(aersn_registry, q == 1 & n == 1000 & prob == 0.99)
  ref1 <- aersn_reference(1, n = 1000, draws = 40000, seed = 20260913,
                          cache = FALSE)
  got1 <- unname(quantile(ref1, 0.99))
  se1 <- sqrt(reg1$mcse^2 + unname(aersn_mcse(ref1, 0.99))^2)
  expect_lt(abs(got1 - reg1$cv), 4 * se1 + 0.005)
})

test_that("p-values, critical values and tables behave as documented", {
  ref <- aersn_reference(1, n = 50, draws = 2000, seed = 3, cache = FALSE)
  cv <- aersn_critical_value(ref, 0.95)
  expect_equal(cv, unname(quantile(ref, 0.95)))
  p <- aersn_pvalue(ref, cv)
  expect_true(p$p.value <= 0.06 && p$p.value >= 0.04)
  expect_equal(p$resolution, 1 / 2000)
  expect_equal(aersn_pvalue(ref, 100)$p.value, 0)
  pg <- aersn_pvalue(ref, 0, alternative = "greater")$p.value
  expect_true(abs(pg - 0.5) < 0.05)
  expect_error(aersn_pvalue(aersn_reference(2, n = 20, draws = 50, seed = 1,
                                            cache = FALSE), 1,
                            alternative = "less"), "q = 1")
  an <- aersn_reference(1, grid = "continuous")
  expect_equal(aersn_pvalue(an, 1.70577645482396)$p.value, 0.05,
               tolerance = 1e-7)
  expect_equal(aersn_pvalue(an, 1.70577645482396, "greater")$p.value, 0.025,
               tolerance = 1e-7)
  expect_equal(aersn_pvalue(an, -1.70577645482396, "less")$p.value, 0.025,
               tolerance = 1e-7)
  expect_true(all(is.na(aersn_mcse(an))))
  ## table reference
  tab <- aersn_table_reference(2, 200, c(0.9, 0.95, 0.99), c(2.2, 2.6, 3.4))
  expect_equal(unname(quantile(tab, 0.95)), 2.6)
  expect_error(quantile(tab, 0.975), "only at")
  expect_equal(aersn_pvalue(tab, 2.0)[c("lower", "upper")],
               list(lower = 0.1, upper = 1))
  expect_equal(aersn_pvalue(tab, 2.7)[c("lower", "upper")],
               list(lower = 0.01, upper = 0.05))
  expect_equal(aersn_pvalue(tab, 5)[c("lower", "upper")],
               list(lower = 0, upper = 0.01))
  expect_output(print(ref), "matched-grid")
  expect_output(print(an), "closed-form")
  expect_output(print(tab), "tabulated")
})

test_that("the bundled registry can be used as a reference", {
  set.seed(14)
  fit <- aersn_mean(matrix(rnorm(400), 200, 2))
  tt <- aersn_test(fit, reference = "registry")
  reg <- subset(aersn_registry, q == 2 & n == 200 & prob == 0.95)
  expect_equal(tt$critical.value, reg$cv)
  expect_error(aersn_test(aersn_mean(rnorm(37)), reference = "registry"),
               "no entry")
})
