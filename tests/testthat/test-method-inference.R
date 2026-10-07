test_that("all methods run at every exercised dimension", {
  for (cfg in list(list(q = 1L, n = 200L), list(q = 2L, n = 200L),
                   list(q = 3L, n = 200L), list(q = 5L, n = 300L),
                   list(q = 10L, n = 400L))) {
    Y <- sim_var1(cfg$n, cfg$q, phi = 0.3, seed = 100 + cfg$q)
    fit <- aersn_mean(Y)
    for (m in aersn_methods()) {
      tt <- aersn_test(fit, method = m, draws = 120, seed = 1)
      expect_true(is.finite(unname(tt$statistic)), info = paste(m, cfg$q))
      expect_true(is.finite(tt$critical.value))
      expect_true(tt$p.value >= 0 && tt$p.value <= 1)
      expect_identical(tt$aersn_method, m)
      ci <- confint(fit, method = m, draws = 120, seed = 1)
      expect_equal(nrow(ci), cfg$q)
      expect_true(all(ci[, "lower"] < ci[, "upper"]))
    }
  }
})

test_that("region membership agrees with the test decision for every method", {
  fit <- aersn_mean(sim_var1(200, 3, phi = 0.4, seed = 110))
  cand <- rbind(c(0, 0, 0), fit$estimate, fit$estimate + c(0.4, -0.3, 0.2))
  for (m in aersn_methods()) {
    reg <- aersn_region(fit, method = m, draws = 200, seed = 2)
    inside <- aersn_contains(reg, cand)
    for (i in seq_len(nrow(cand))) {
      tt <- aersn_test(fit, null = cand[i, ], method = m, draws = 200, seed = 2)
      expect_identical(unname(inside[i]), !tt$reject, info = paste(m, i))
    }
    ## The centre is always inside, with gauge zero.
    expect_equal(unname(attr(inside, "gauge")[2L]), 0, tolerance = 1e-12)
    ## Boundary points have region gauge one.
    B <- aersn_slice(reg, c(1L, 2L), n_angles = 12L)
    full <- cbind(B, fit$estimate[3L])
    expect_equal(unname(attr(aersn_contains(reg, full), "gauge")),
                 rep(1, 12), tolerance = 1e-7, info = m)
  }
})

test_that("the support function equals the simultaneous interval endpoint", {
  fit <- aersn_mean(sim_var1(200, 3, phi = 0.4, seed = 111))
  A <- rbind(c(1, 0, 0), c(1, -1, 0), c(0.3, 0.2, -0.5))
  for (m in aersn_methods()) {
    reg <- aersn_region(fit, method = m, draws = 200, seed = 3)
    ct <- aersn_contrast(fit, A, method = m, draws = 200, seed = 3)
    expect_equal(ct$upper, aersn_support(reg, A), tolerance = 1e-9, info = m)
    expect_equal(ct$lower, -aersn_support(reg, -A), tolerance = 1e-9)
    ci <- confint(fit, method = m, draws = 200, seed = 3)
    expect_equal(unclass(ci)[, ], unclass(confint(reg))[, ], tolerance = 1e-9)
  }
})

test_that("interval types use the documented reference dimension and ordering", {
  fit <- aersn_mean(sim_var1(300, 3, phi = 0.3, seed = 112))
  for (m in aersn_methods()) {
    cs <- confint(fit, method = m, draws = 400, seed = 4)
    cj <- confint(fit, parm = 1:2, type = "joint", method = m, draws = 400,
                  seed = 4)
    cmg <- confint(fit, type = "marginal", method = m, draws = 400, seed = 4)
    expect_identical(attr(cs, "dimension"), 3L)
    expect_identical(attr(cj, "dimension"), 2L)
    expect_identical(attr(cmg, "dimension"), 1L)
    w <- function(x) unclass(x)[, "upper"] - unclass(x)[, "lower"]
    ## Wider reference dimension, wider interval.
    expect_true(all(w(cs)[1:2] > w(cj)), info = m)
    expect_true(all(w(cj) > w(cmg)[1:2]), info = m)
    expect_output(print(cmg), "do not provide simultaneous")
  }
})

test_that("joint inference on a target rebuilds the method on that target", {
  fit <- aersn_mean(sim_var1(250, 3, phi = 0.5, seed = 113))
  A <- rbind(c(1, -1, 0))
  ## For the hull and for every normalizer that is a linear functional of
  ## outer products, rebuilding on the target agrees with transforming the
  ## full-dimensional normalizer.
  for (m in c("shao", "fixedb", "ewc")) {
    nz_full <- aersn_normalizer(fit, m)
    tgt <- aersn_target(fit, A, names = "d")
    nz_t <- aersn_normalizer(tgt, m)
    expect_equal(unname(nz_t$matrix), unname(A %*% nz_full$matrix %*% t(A)),
                 tolerance = 1e-9, info = m)
  }
  ## HAC with a fixed bandwidth also transforms exactly.
  nz_full <- aersn_hac_lrv(fit, bandwidth = 6)
  nz_t <- aersn_hac_lrv(aersn_target(fit, A, names = "d"), bandwidth = 6)
  expect_equal(unname(nz_t$matrix), unname(A %*% nz_full$matrix %*% t(A)),
               tolerance = 1e-9)
  ## The LDL factorization is recomputed for the target and does not.
  ldl_full <- aersn_ldl_normalizer(fit)
  ldl_t <- aersn_ldl_normalizer(aersn_target(fit, A, names = "d"))
  expect_false(isTRUE(all.equal(unname(ldl_t$matrix),
                                unname(A %*% ldl_full$matrix %*% t(A)))))
  ## Hull coordinate projections are unchanged by selection, as before.
  sel <- rbind(c(1, 0, 0), c(0, 1, 0))
  expect_equal(unname(aersn_target(fit, sel)$hull$coordinate_ranges),
               unname(fit$hull$coordinate_ranges[1:2]))
})

test_that("statistics are invariant under nonsingular reparameterization where the method guarantees it", {
  Y <- sim_var1(250, 3, phi = 0.4, seed = 114)
  fit <- aersn_mean(Y)
  set.seed(115)
  Q <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  D <- diag(c(0.01, 1, 100))
  S <- diag(3); S[1, 2] <- 0.8; S[3, 1] <- -0.4
  z0 <- matrix(sqrt(fit$n) * fit$estimate, nrow = 1L)
  for (H in list(rotation = Q, rescaling = D, shear = S)) {
    fitH <- aersn_mean(Y %*% t(H))
    zH <- matrix(sqrt(fitH$n) * fitH$estimate, nrow = 1L)
    ## Invariant: the hull, and every Wald statistic whose normalizer is a
    ## linear functional of outer products with tuning that does not depend
    ## on the coordinate system.
    for (m in c("hull", "shao", "fixedb", "ewc")) {
      expect_equal(
        normalizer_statistic(aersn_normalizer(fitH, m), zH),
        normalizer_statistic(aersn_normalizer(fit, m), z0),
        tolerance = 1e-7, info = m)
    }
    expect_equal(
      normalizer_statistic(aersn_hac_lrv(fitH, bandwidth = 6), zH),
      normalizer_statistic(aersn_hac_lrv(fit, bandwidth = 6), z0),
      tolerance = 1e-7)
  }
  ## Not invariant, and recorded as such: the LDL factorization depends on the
  ## coordinate system, and an automatic bandwidth is selected from the
  ## transformed contributions.
  fitQ <- aersn_mean(Y %*% t(Q))
  zQ <- matrix(sqrt(fitQ$n) * fitQ$estimate, nrow = 1L)
  expect_false(isTRUE(all.equal(
    normalizer_statistic(aersn_ldl_normalizer(fitQ), zQ),
    normalizer_statistic(aersn_ldl_normalizer(fit), z0))))
  expect_false(isTRUE(all.equal(
    aersn_hac_lrv(fitQ, bandwidth = "andrews")$tuning$bandwidth,
    aersn_hac_lrv(fit, bandwidth = "andrews")$tuning$bandwidth)))
  ## A diagonal rescaling does leave the componentwise statistic unchanged,
  ## because it does not mix coordinates.
  fitD <- aersn_mean(Y %*% t(D))
  zD <- matrix(sqrt(fitD$n) * fitD$estimate, nrow = 1L)
  expect_equal(normalizer_statistic(aersn_ldl_normalizer(fitD), zD),
               normalizer_statistic(aersn_ldl_normalizer(fit), z0),
               tolerance = 1e-6)
})

test_that("one-sided alternatives are restricted to the hull method", {
  fit <- aersn_mean(as.numeric(sim_var1(200, 1, seed = 116)))
  expect_s3_class(aersn_test(fit, 0, alternative = "greater",
                             reference = "continuous"), "htest")
  for (m in setdiff(aersn_methods(), "hull")) {
    expect_error(aersn_test(fit, 0, alternative = "greater", method = m),
                 "method")
  }
})

test_that("two-dimensional region geometry is method specific", {
  fit <- aersn_mean(sim_var1(200, 2, phi = 0.4, seed = 117))
  hull <- aersn_region(fit, method = "hull", draws = 200, seed = 5)
  ell <- aersn_region(fit, method = "shao", draws = 200, seed = 5)
  expect_identical(hull$family, "geometric")
  expect_identical(ell$family, "quadratic")
  V <- aersn_vertices(hull)
  expect_true(nrow(V) >= 3)
  expect_error(aersn_vertices(ell), "ellipsoid")
  ## The ellipse area formula agrees with the polygon area of its boundary.
  B <- aersn_projection(ell, c(1L, 2L))
  expect_equal(region_area_2d(ell), polygon_area(B), tolerance = 1e-3)
  expect_equal(region_area_2d(hull), polygon_area(V), tolerance = 1e-10)
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_invisible(plot(ell))
  expect_invisible(plot(hull, null = c(0, 0)))
})
