test_that("the cosine basis is orthonormal and orthogonal to the constant", {
  for (n in c(40L, 101L, 250L)) {
    nu <- min(12L, n - 1L)
    Phi <- ewc_basis(n, nu)
    expect_equal(dim(Phi), c(n, nu))
    expect_equal(crossprod(Phi), diag(nu), tolerance = 1e-10)
    expect_lt(max(abs(colSums(Phi))), 1e-10)
  }
})

test_that("the EWC matrix matches the replication implementation", {
  for (q in c(1L, 2L, 3L, 5L)) {
    Y <- sim_var1(200, q, phi = 0.4, seed = 80 + q)
    fit <- aersn_mean(Y)
    for (nu in c(q, 10L, 24L)) {
      nz <- aersn_ewc_lrv(fit, nu = nu)
      expect_equal(unname(nz$matrix), unname(ref_ewc_matrix(fit$psi, nu)),
                   tolerance = 1e-12, info = paste(q, nu))
      ## The replication code passes the raw series; the basis is orthogonal
      ## to the constant, so demeaning cannot change the estimate.
      expect_equal(unname(nz$matrix), unname(ref_ewc_matrix(Y, nu)),
                   tolerance = 1e-10)
      expect_equal(aersn_ewc_lrv(fit, nu = nu, center = FALSE)$matrix,
                   nz$matrix, tolerance = 1e-10)
    }
  }
})

test_that("the number of cosine terms follows the documented rule and limits", {
  for (n in c(120L, 300L, 1000L)) {
    fit <- aersn_mean(sim_var1(n, 2, seed = 85))
    nz <- aersn_ewc_lrv(fit)
    expect_identical(nz$tuning$nu, as.integer(floor(0.4 * n^(2 / 3))))
    expect_identical(nz$tuning$rule, "floor(0.4 n^(2/3))")
    expect_equal(nz$diagnostics$smoothing_fraction, nz$tuning$nu / n)
  }
  fit <- aersn_mean(sim_var1(60, 3, seed = 86))
  expect_error(aersn_ewc_lrv(fit, nu = 2), "at least as many cosine terms")
  expect_s3_class(aersn_ewc_lrv(fit, nu = 3), "aersn_normalizer")
  expect_error(aersn_ewc_lrv(fit, nu = 60), "at most n - 1")
  expect_s3_class(aersn_ewc_lrv(fit, nu = 59), "aersn_normalizer")
  expect_error(aersn_ewc_lrv(fit, nu = 4.5), "whole number")
})

test_that("the EWC reference law is the scaled F law of the replication code", {
  for (q in c(1L, 2L, 3L, 5L)) {
    for (nu in c(q, q + 3L, 30L)) {
      ref <- aersn_parametric_reference("ewc_F", df = q, nu = nu)
      expect_equal(aersn_critical_value(ref, 0.95),
                   ref_ewc_critical(q, nu, 0.95), tolerance = 1e-10)
      expect_equal(aersn_critical_value(ref, 0.99),
                   ref_ewc_critical(q, nu, 0.99), tolerance = 1e-10)
      ## The p-value inverts the same transformation.
      stat <- ref_ewc_critical(q, nu, 0.95)
      expect_equal(aersn_pvalue(ref, stat)$p.value, 0.05, tolerance = 1e-8)
    }
  }
  expect_error(aersn_parametric_reference("ewc_F", df = 3, nu = 2), "nu >= df")
  expect_error(aersn_parametric_reference("ewc_F", df = 2), "`nu` is required")
})

test_that("a linear restriction uses the dimension of the restriction", {
  fit <- aersn_mean(sim_var1(300, 3, seed = 87))
  nu <- aersn_ewc_lrv(fit)$tuning$nu
  A <- rbind(c(1, -1, 0))
  tt <- aersn_test(fit, null = 0, contrast = A, method = "ewc")
  expect_equal(tt$critical.value, ref_ewc_critical(1, nu, 0.95),
               tolerance = 1e-10)
  A2 <- rbind(c(1, -1, 0), c(0, 1, -1))
  tt2 <- aersn_test(fit, null = c(0, 0), contrast = A2, method = "ewc")
  expect_equal(tt2$critical.value, ref_ewc_critical(2, nu, 0.95),
               tolerance = 1e-10)
  ## The full joint test uses q = 3.
  tt3 <- aersn_test(fit, method = "ewc")
  expect_equal(tt3$critical.value, ref_ewc_critical(3, nu, 0.95),
               tolerance = 1e-10)
  ## Restricting to a contrast is not the same as projecting the joint region.
  expect_true(tt$critical.value < tt3$critical.value)
})

test_that("the chi-squared reference law is used for HAC only", {
  ref <- aersn_parametric_reference("chisq", df = 4)
  expect_equal(aersn_critical_value(ref, 0.95), stats::qchisq(0.95, 4))
  expect_equal(aersn_pvalue(ref, 9.49)$p.value,
               stats::pchisq(9.49, 4, lower.tail = FALSE))
  expect_true(all(is.na(aersn_mcse(ref))))
  expect_error(aersn_pvalue(ref, 1, alternative = "greater"),
               "One-sided alternatives are not defined")
  expect_output(print(ref), "chi-squared")
})
