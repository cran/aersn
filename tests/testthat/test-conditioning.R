test_that("QR preconditioning preserves gauges and maps dual directions back", {
  for (q in c(2L, 3L, 5L)) {
    G <- random_path(60, q, seed = 701 + q)$G
    z <- seq_len(q) / q
    value <- aersn:::gauge_lp(G, z)$value
    for (epsilon in c(1e-6, 1e-8, 1e-10)) {
      H <- diag(q)
      H[q, ] <- c(rep(1, q - 1), epsilon)
      transformed <- G %*% t(H)
      error <- as.numeric(H %*% z)
      result <- aersn:::gauge_lp(transformed, error)
      expect_equal(result$value, value, tolerance = 1e-4)
      expect_equal(sum(result$direction * error), value, tolerance = 1e-4)
      projection <- as.numeric(transformed %*% result$direction)
      expect_equal(diff(range(projection)), 1, tolerance = 1e-4)
      expect_gte(min(projection) - result$offset, -1e-4)
      expect_lte(max(projection) - result$offset, 1 + 1e-4)
    }
    # Units and permutations exercise both the range scale and QR pivot.
    H <- diag(10^seq(-9, 9, length.out = q))[q:1, , drop = FALSE]
    result <- aersn:::gauge_lp(G %*% t(H), as.numeric(H %*% z))
    expect_equal(result$value, value, tolerance = 1e-8)
    expect_equal(sum(result$direction * as.numeric(H %*% z)),
                 value, tolerance = 1e-8)
  }
})

test_that("ill-conditioned joint tests agree with simultaneous region projections", {
  set.seed(3)
  n <- 60
  psi0 <- matrix(rnorm(2 * n), n, 2)
  H <- rbind(c(1, 0), c(1, 1e-11))
  psi <- psi0 %*% t(H)
  fit <- aersn(colMeans(psi), psi)
  null <- fit$estimate - as.numeric(H %*% c(0.3, 4)) / sqrt(n)
  base <- aersn(colMeans(psi0), psi0)
  expected <- aersn_gauge(base, base$estimate - c(0.3, 4) / sqrt(n))
  ref <- aersn_reference(fit, draws = 200, seed = 1)
  result <- aersn_test(fit, null = null, reference = ref)
  expect_equal(unname(result$statistic), expected, tolerance = 1e-3)
  expect_true(result$reject)
  region <- aersn_region(fit, reference = ref)
  expect_false(as.logical(aersn_contains(region, null)))
  ci <- aersn_contrast(fit, c(-1, 1), reference = ref)
  contrast_null <- sum(c(-1, 1) * null)
  expect_true(contrast_null < ci$lower || contrast_null > ci$upper)
})

test_that("QR handles zero errors, repeated nodes and rejects singular paths", {
  G <- random_path(20, 3, seed = 731)$G
  duplicated <- G[rep(seq_len(nrow(G)), each = 2), ]
  expect_equal(aersn:::gauge_lp(duplicated, c(1, 2, 3))$value,
               aersn:::gauge_lp(G, c(1, 2, 3))$value, tolerance = 1e-8)
  expect_equal(aersn:::gauge_lp(G, rep(0, 3))$value, 0)
  bad <- cbind(G[, 1], G[, 1])
  expect_error(aersn:::gauge_lp(bad, c(1, 1)), "rank deficient|ill-conditioned")
})
