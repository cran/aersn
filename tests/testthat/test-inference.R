test_that("aersn() validates its inputs", {
  psi <- matrix(rnorm(40), 20, 2)
  expect_error(aersn(c(0, 0, 0), psi), "one column per parameter")
  expect_error(aersn(c(0, 0), psi, n = 19), "one row per observation")
  expect_error(aersn(c(NA, 0), psi), "finite")
  fit <- aersn(c(0.1, -0.2), psi, names = c("a", "b"))
  expect_equal(coef(fit), c(a = 0.1, b = -0.2))
  expect_equal(nobs(fit), 20L)
  expect_output(print(fit), "increment hull")
  expect_error(aersn_test(fit, null = c(0, NA)), "finite")
  expect_error(aersn_test(fit, null = 0), "length 2")
  expect_error(aersn_test(fit, level = 1), "strictly between")
})

test_that("test statistics, regions and intervals are mutually consistent", {
  set.seed(15)
  n <- 120
  Y <- matrix(rnorm(n * 2), n, 2)
  fit <- aersn_mean(Y)
  ref <- aersn_reference(2, n = n, draws = 400, seed = 2)
  tt <- aersn_test(fit, null = c(0, 0), reference = ref)
  expect_equal(unname(tt$statistic), aersn_gauge(fit, c(0, 0)))
  expect_equal(unname(tt$statistic),
               aersn:::gauge_lp(fit$path$G, sqrt(n) * fit$estimate)$value)
  reg <- aersn_region(fit, reference = ref)
  cv <- tt$critical.value
  ## membership equals non-rejection
  cand <- rbind(c(0, 0), fit$estimate, fit$estimate + c(0.5, -0.5))
  inside <- aersn_contains(reg, cand)
  stats <- aersn_gauge(fit, cand)
  expect_equal(as.logical(inside), stats <= cv)
  expect_equal(attr(inside, "gauge"), stats / cv)
  ## vertices lie on the boundary; the center is inside
  V <- aersn_vertices(reg)
  expect_equal(aersn_gauge(reg, V), rep(1, nrow(V)), tolerance = 1e-8)
  expect_equal(aersn_gauge(reg, fit$estimate), 0)
  ## support function of the region equals the simultaneous contrast bound
  A <- rbind(c(1, 0), c(0, 1), c(1, -1))
  ct <- aersn_contrast(fit, A, reference = ref)
  expect_equal(ct$upper, aersn_support(reg, A))
  expect_equal(ct$lower, -aersn_support(reg, -A))
  ci <- confint(fit, reference = ref)
  expect_equal(unname(ci[, "upper"]), ct$upper[1:2])
  expect_equal(unclass(confint(reg))[, ], unclass(ci)[, ])
  ## the boundary point achieving the projection satisfies T = c
  j <- 1
  proj <- fit$path$G[, j]
  xstar <- fit$path$G[which.max(proj), ] - fit$path$G[which.min(proj), ]
  vstar <- fit$estimate + reg$scale * xstar
  expect_equal(aersn_gauge(fit, vstar), cv, tolerance = 1e-8)
  expect_equal(unname(vstar[j]), unname(ci[j, "upper"]))
  ## polygon area is positive and region printing works
  expect_output(print(reg), "area")
  expect_output(print(tt), "critical value")
  expect_output(print(ci), "Simultaneous")
  expect_output(print(ct), "contrasts")
})

test_that("interval types are distinguished and ordered", {
  set.seed(16)
  fit <- aersn_mean(matrix(rnorm(300), 100, 3))
  cs <- confint(fit, draws = 1500, seed = 4)
  cj <- confint(fit, parm = 1:2, type = "joint", draws = 1500, seed = 4)
  cm <- confint(fit, type = "marginal", draws = 6000, seed = 4)
  expect_equal(attr(cs, "dimension"), 3L)
  expect_equal(attr(cj, "dimension"), 2L)
  expect_equal(attr(cm, "dimension"), 1L)
  expect_gt(attr(cs, "critical.value"), attr(cj, "critical.value"))
  expect_gt(attr(cj, "critical.value"), attr(cm, "critical.value"))
  w <- function(ci) unclass(ci)[, "upper"] - unclass(ci)[, "lower"]
  expect_true(all(w(cs)[1:2] > w(cj)))
  expect_true(all(w(cj) > w(cm)[1:2]))
  expect_output(print(cm), "do not provide simultaneous")
  expect_error(confint(fit, parm = "zz"), "Unknown parameter")
  expect_error(confint(fit, parm = c(1, 1)), "distinct")
  ct <- aersn_contrast(fit, c(1, -1, 0), type = "marginal", draws = 6000,
                       seed = 4)
  expect_equal(attr(ct, "dimension"), 1L)
  expect_error(aersn_contrast(fit, rbind(c(1, 0, 0), c(2, 0, 0)),
                              type = "joint", draws = 100), "full row rank")
})

test_that("linear restrictions are tested as point nulls of the transformed target", {
  set.seed(17)
  fit <- aersn_mean(matrix(rnorm(300), 100, 3), names = c("a", "b", "c"))
  A <- rbind(c(1, -1, 0), c(0, 1, -1))
  ref <- aersn_reference(2, n = 100, draws = 300, seed = 5)
  t1 <- aersn_test(fit, null = c(0, 0), contrast = A, reference = ref)
  t2 <- aersn_test(aersn_target(fit, A), null = c(0, 0), reference = ref)
  expect_equal(t1$statistic, t2$statistic)
  expect_equal(t1$p.value, t2$p.value)
  ## function target with Jacobian
  h <- function(b) c(b[1] - b[2], b[2] - b[3])
  t3 <- aersn_test(aersn_target(fit, h, jacobian = A), reference = ref)
  expect_equal(t3$statistic, t1$statistic)
  expect_error(aersn_target(fit, rbind(c(1, 0, 0), c(1, 0, 0))),
               "full row rank")
  expect_error(aersn_target(fit, h), "jacobian")
  expect_error(aersn_target(fit, matrix(1, 1, 2)), "columns")
})

test_that("scalar one-sided tests use the signed ratio", {
  set.seed(18)
  x <- rnorm(80, mean = 0.3)
  fit <- aersn_mean(x)
  two <- aersn_test(fit, 0, reference = "continuous")
  gr <- aersn_test(fit, 0, alternative = "greater", reference = "continuous")
  le <- aersn_test(fit, 0, alternative = "less", reference = "continuous")
  expect_equal(unname(gr$statistic), unname(two$statistic))
  expect_equal(gr$p.value, two$p.value / 2)
  expect_equal(le$p.value, 1 - two$p.value / 2)
  expect_equal(gr$critical.value, aersn_scalar_quantile(0.90))
  expect_equal(le$critical.value, -aersn_scalar_quantile(0.90))
  expect_true(le$critical.value < 0)
  fit2 <- aersn_mean(matrix(rnorm(80), 40, 2))
  expect_error(aersn_test(fit2, alternative = "greater"), "q = 1")
})

test_that("projections and slices of higher-dimensional regions are consistent", {
  set.seed(19)
  fit <- aersn_mean(matrix(rnorm(450), 150, 3))
  ref <- aersn_reference(3, n = 150, draws = 300, seed = 6)
  reg <- aersn_region(fit, reference = ref)
  P <- aersn_projection(reg, c(1, 3))
  ## the projection is the region of the two-dimensional target with c_3
  sub <- aersn_target(fit, diag(3)[c(1, 3), ])
  V <- aersn:::hull_vertices_2d(sub$path$G)
  P2 <- sweep(V * reg$scale, 2, fit$estimate[c(1, 3)], `+`)
  expect_equal(P, P2, ignore_attr = TRUE)
  ## slice points are on the boundary of the region
  S <- aersn_slice(reg, c(1, 2), n_angles = 24)
  full <- cbind(S, fit$estimate[3])
  expect_equal(aersn_gauge(reg, full), rep(1, 24), tolerance = 1e-8)
  expect_error(aersn_vertices(reg), "q = 2")
  expect_error(aersn_projection(reg, c(1, 1)), "distinct")
  pdf(NULL); on.exit(dev.off())
  expect_invisible(plot(reg))
  expect_invisible(plot(reg, type = "slice", null = c(0, 0, 0)))
  expect_invisible(plot(fit$path))
  reg2 <- aersn_region(aersn_mean(matrix(rnorm(100), 50, 2)), draws = 100,
                       seed = 1)
  expect_invisible(plot(reg2, null = c(1, 1)))
  reg1 <- aersn_region(aersn_mean(rnorm(50)), draws = 100, seed = 1)
  expect_invisible(plot(reg1, null = 0))
  s <- summary(fit, reference = ref)
  expect_output(print(s), "Simultaneous")
})
