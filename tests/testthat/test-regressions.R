test_that("reference caching respects uncertainty settings and profile nodes", {
  aersn_clear_cache()
  a <- aersn_reference(1, n = 40, draws = 400, seed = 928, batches = 4)
  b <- aersn_reference(1, n = 40, draws = 400, seed = 928, batches = 20)
  c <- aersn_reference(1, n = 40, draws = 400, seed = 928, batches = 20,
                        cache = FALSE)
  expect_identical(b$batches, 20L)
  expect_identical(a$draws, b$draws)
  expect_identical(aersn_mcse(b, .95), aersn_mcse(c, .95))
  expect_false(identical(aersn_mcse(a, .95), aersn_mcse(b, .95)))
  expect_true(all(is.na(aersn_mcse(aersn_reference(1, n = 10, draws = 10)))))
})

test_that("reference grids require enough positive increments, allowing flat segments", {
  expect_error(aersn_reference(1, nodes = c(rep(0, 20), 1), draws = 20),
               "degenerate bridge")
  expect_error(aersn_reference(2, nodes = c(0, 0, .5, 1, 1), draws = 20),
               "positive grid increments")
  nodes <- c(0, 0, .2, .6, 1, 1)
  ref <- aersn_reference(2, nodes = nodes, draws = 30, cache = FALSE)
  expect_length(ref$draws, 30)
  expect_true(all(is.finite(ref$draws)))
  expect_identical(ref$lp_failures, 0L)
  expect_error(aersn_reference(1, n = 30, nodes = nodes), "length\\(nodes\\)")
  expect_error(aersn_reference(1, n = 30, draws = 2^31), "integer range")
  expect_error(aersn_reference(1, n = 30, cores = 3), "cores")
  expect_error(aersn_reference(1, n = 30, cache = NA), "TRUE or FALSE")
})

test_that("failed simulated draws are not discarded or cached, and RNG is restored", {
  aersn_clear_cache()
  set.seed(100)
  before <- .Random.seed
  count <- 0L
  original <- gauge_lp
  local_mocked_bindings(gauge_lp = function(...) {
    count <<- count + 1L
    if (count == 3L) stop("injected numerical failure")
    original(...)
  })
  expect_error(aersn_reference(2, n = 10, draws = 5),
               "1 of 5 reference draws failed numerically")
  expect_identical(.Random.seed, before)
  expect_identical(aersn_clear_cache(), 0L)
})

test_that("OLS na.exclude and na.omit use exactly the same fitted rows", {
  set.seed(9)
  d <- data.frame(x = rnorm(100), y = rnorm(100))
  d$x[21] <- NA
  expect_warning(a <- aersn_lm(lm(y ~ x, d, na.action = na.exclude)),
                 "treated as consecutive")
  expect_warning(b <- aersn_lm(lm(y ~ x, d, na.action = na.omit)),
                 "treated as consecutive")
  expect_equal(a$estimate, b$estimate)
  expect_equal(a$psi, b$psi)
  expect_equal(a$path$G, b$path$G)
  expect_identical(a$n, 99L)
  expect_error(aersn_lm(glm(y ~ x, data = d)), "single-response")
  d$x[21] <- 0
  expect_error(aersn_lm(lm(cbind(y, x) ~ 1, d)), "single-response")
})

test_that("invalid coordinate selectors do not silently change the target", {
  set.seed(5)
  fit <- aersn_mean(matrix(rnorm(60), 30, 2), names = c("a", "b"))
  for (idx in list(1.5, NA_real_, Inf, numeric(), c(1, 1))) {
    expect_error(confint(fit, parm = idx, draws = 20), "integer coordinates")
    expect_error(resolve_target(idx, c(a = 1, b = 2), 2), "integer coordinates")
  }
  expect_error(confint(fit, parm = c("a", "a")), "integer coordinates")
  expect_error(resolve_target(matrix(1, 1, 3), c(1, 2), 2), "columns")
  expect_error(aersn_target(fit, matrix(numeric(), 0, 2)), "at least one row")
  expect_error(aersn_path(matrix(numeric(), 5, 0)), "one column")
  expect_error(aersn_path(1), "two observations")
  expect_error(aersn_profile(NULL, 5, tol = NA_real_), "tol")
  reg <- aersn_region(fit, draws = 20)
  expect_error(aersn_projection(reg, c(1.5, 2)), "integer coordinates")
  expect_error(aersn_slice(reg, n_angles = 2), "n_angles")
  expect_error(aersn_contains(reg, c(0, 0), tol = NA_real_), "tol")
})

test_that("joint intervals for an identified subset allow a redundant full vector", {
  set.seed(10)
  z <- rnorm(40)
  expect_warning(fit <- aersn_mean(cbind(a = z, b = z)), "rank deficient")
  ref <- aersn_reference(1, n = 40, draws = 100)
  ci <- confint(fit, parm = 1, type = "joint", reference = ref)
  direct <- confint(aersn_mean(z), reference = ref)
  expect_equal(as.numeric(ci), as.numeric(direct))
  expect_error(confint(fit, type = "simultaneous", reference = ref),
               "full dimensional")
})

test_that("tabulated references validate quantiles, probabilities and standard errors", {
  for (values in list(NA_real_, Inf, -1, 0)) {
    expect_error(aersn_table_reference(1, 10, .95, values), "finite positive")
  }
  expect_error(aersn_table_reference(2, 2, .95, 2), "n >= q")
  expect_error(aersn_table_reference(1, 10, c(.95, .95), c(1, 2)), "distinct")
  expect_error(aersn_table_reference(1, 10, c(.9, .95), c(2, 1)), "nondecreasing")
  for (se in list(-1, Inf, c(.1, .2))) {
    expect_error(aersn_table_reference(1, 10, .95, 2, mcse = se), "mcse")
  }
  ref <- aersn_table_reference(1, 10, c(.99, .95), c(3, 2), c(.2, .1))
  expect_equal(unname(aersn_mcse(ref, .95)), .1)
  expect_equal(aersn_pvalue(ref, -2.5), aersn_pvalue(ref, 2.5))
  expect_error(aersn_mcse(ref, 1), "strictly between")
})

test_that("one-sided probabilities and critical values use the same symmetric law", {
  ref <- aersn_reference(1, n = 40, draws = 400, seed = 17)
  for (z in c(-2, -.5, 0, .5, 2)) {
    p2 <- aersn_pvalue(ref, z)$p.value
    pg <- aersn_pvalue(ref, z, "greater")$p.value
    pl <- aersn_pvalue(ref, z, "less")$p.value
    expect_equal(pg + pl, 1)
    expect_equal(pg, aersn_pvalue(ref, -z, "less")$p.value)
    expect_equal(min(pg, pl), p2 / 2)
  }
  ## These are the rejection boundaries implied by the one-sided tails.
  cv <- aersn_critical_value(ref, .9)
  expect_lt(abs(aersn_pvalue(ref, cv, "greater")$p.value - .05), 1 / 400)
  set.seed(8)
  fit <- aersn_mean(rnorm(40))
  gt <- aersn_test(fit, alternative = "greater", reference = ref)
  lt <- aersn_test(fit, alternative = "less", reference = ref)
  expect_equal(gt$critical.value, cv)
  expect_equal(lt$critical.value, -cv)
  zero <- aersn_pvalue(ref, max(ref$draws) + 1)
  expect_equal(zero$p.value, 0)
  expect_equal(zero$conf.int[2], 1 - .025^(1 / 400))
  expect_gt(zero$conf.int[2], zero$resolution)
  fit$estimate[] <- 100
  expect_output(print(aersn_test(fit, reference = ref)), "95% Monte Carlo interval")
})

test_that("analytic tails remain accurate beyond CDF subtraction precision", {
  ref <- aersn_reference(1, grid = "continuous")
  p <- aersn_pvalue(ref, 20)$p.value
  expect_gt(p, 0)
  expect_equal(p, 2 * pi * 20^2 * besselK(pi * 20, 1), tolerance = 1e-12)
  expect_equal(aersn_scalar_cdf(1e200), 1)
  expect_equal(aersn_scalar_density(1e200), 0)
  expect_error(aersn_scalar_cdf(1, tol = 0), "tol")
  expect_error(aersn_scalar_density(1, tol = NA_real_), "tol")
})

test_that("invalid weighting and Hessian matrices fail before constructing scores", {
  set.seed(2)
  s <- matrix(rnorm(80), 40, 2)
  s <- sweep(s, 2, colMeans(s))
  expect_error(aersn_gmm(s, diag(2), c(0, 0), weight = diag(c(1, Inf))), "finite")
  expect_error(aersn_mle(s, matrix(c(1, 0, 1, 1), 2), c(0, 0)), "symmetric")
  expect_error(aersn_mle(s, -diag(2), c(0, 0)), "positive definite")
})

test_that("changing coordinate units does not remove identifiable parameters", {
  set.seed(18)
  y <- matrix(rnorm(300), 100, 3)
  fit <- aersn_mean(y)
  units <- c(1e-9, 1, 1e9)
  scaled <- aersn_mean(sweep(y, 2, units, `*`))
  expect_true(scaled$hull$full_rank)
  expect_lt(scaled$hull$raw_rank, scaled$q)
  expect_equal(aersn_gauge(scaled, c(0, 0, 0)),
               aersn_gauge(fit, c(0, 0, 0)), tolerance = 1e-8)
  ## A true linear dependence must still be detected after rescaling.
  redundant <- cbind(y[, 1:2], y[, 1] + y[, 2])
  expect_warning(deg <- aersn_mean(sweep(redundant, 2, units, `*`)),
                 "rank deficient")
  expect_equal(deg$hull$rank, 2)
})

test_that("excessive analytic series are rejected rather than silently truncated", {
  expect_error(aersn_scalar_cdf(1e3, method = "series"), "too many terms")
  expect_error(aersn_scalar_cdf(1e-10, method = "bessel"), "too many terms")
  expect_error(aersn_scalar_density(1e5, method = "series"), "too many terms")
  expect_equal(aersn_scalar_cdf(1e-10), 1e-10, tolerance = 1e-12)
})

test_that("old or damaged reference draws cannot silently enter inference", {
  ref <- aersn_reference(1, n = 20, draws = 40)
  ref$draws <- ref$draws[-1]
  expect_error(aersn_critical_value(ref), "incomplete")
  ref <- aersn_reference(1, n = 20, draws = 40)
  ref$draws[1] <- Inf
  expect_error(aersn_pvalue(ref, 1), "incomplete")
  ref <- aersn_reference(1, n = 20, draws = 40)
  ref$lp_failures <- 1L
  expect_error(quantile(ref), "failed draws")
})

test_that("region plots accept user labels and restore the graphics layout", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)
  set.seed(12)
  for (q in 1:3) {
    fit <- aersn_mean(matrix(rnorm(40 * q), 40, q))
    reg <- aersn_region(fit, draws = 30)
    old <- graphics::par(c("mfrow", "mar"))
    expect_no_condition(plot(reg, main = "My region", xlab = "First", ylab = "Second"))
    expect_equal(graphics::par(c("mfrow", "mar")), old)
  }
})

test_that("the score profile is invariant to coordinate units and nonsingular maps", {
  set.seed(84)
  s <- matrix(rnorm(240), 80, 3)
  s <- sweep(s, 2, colMeans(s))
  base <- aersn_opg_profile(s)
  A <- matrix(c(1, 2, 0, 0, 1, 1, 0, 0, 1), 3)
  for (ss in list(sweep(s, 2, c(1e-9, 1, 1e9), `*`), s %*% t(A))) {
    mapped <- aersn_opg_profile(ss)
    expect_equal(as.numeric(mapped), as.numeric(base), tolerance = 1e-10)
    expect_equal(attr(mapped, "proportionality"), attr(base, "proportionality"),
                 tolerance = 1e-10)
  }
})
