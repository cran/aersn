test_that("hull diagnostics detect rank deficiency and degenerate inputs", {
  set.seed(9)
  psi <- matrix(rnorm(100), 50, 2)
  psi3 <- cbind(psi, psi[, 1] + psi[, 2])      # exactly dependent column
  expect_warning(fit <- aersn(c(0, 0, 0), psi3), "rank deficient")
  expect_false(fit$hull$full_rank)
  expect_equal(fit$hull$rank, 2L)
  expect_error(aersn_gauge(fit, c(0, 0, 0)), "not full dimensional")
  expect_error(aersn_test(fit, c(0, 0, 0)), "not full dimensional")
  expect_error(confint(fit, draws = 100), "not full dimensional")
  expect_error(aersn_region(fit, draws = 100), "not full dimensional")
  ## marginal intervals remain available
  ci <- confint(fit, type = "marginal", draws = 500)
  expect_equal(dim(ci), c(3L, 2L))
  ## constant series
  expect_error(aersn_mean(rep(1, 20)), "zero adjusted range")
  ## too few observations
  expect_error(aersn_hull(aersn_path(matrix(rnorm(6), 2, 3))), "q \\+ 1")
})

test_that("hull diagnostics are computed and printed", {
  path <- random_path(60, 3, seed = 10)
  hull <- aersn_hull(path, n_directions = 50)
  expect_true(hull$full_rank)
  expect_equal(hull$rank, 3L)
  expect_equal(hull$n_directions, 53L)
  expect_true(hull$min_projected_range > 0)
  expect_gte(hull$projected_range_spread, 1)
  expect_equal(unname(hull$coordinate_ranges),
               unname(apply(path$G, 2, function(x) diff(range(x)))))
  expect_output(print(hull), "full dimensional")
  expect_output(print(path), "Centered influence path")
})
