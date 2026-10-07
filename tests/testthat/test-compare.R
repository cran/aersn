test_that("the comparison table has the documented columns and order", {
  fit <- aersn_mean(sim_var1(250, 2, phi = 0.4, seed = 140))
  cmp <- aersn_compare(fit, null = c(0, 0), draws = 150, seed = 1)
  expect_s3_class(cmp, "aersn_comparison")
  expect_identical(cmp$method,
                   c("hull", "ldl", "shao", "hac", "fixedb", "ewc"))
  expect_true(all(c("method", "label", "statistic", "statistic_scale",
                    "critical_value", "p_value", "p_mcse", "reject",
                    "half_width", "reference_family", "reference", "tuning",
                    "status") %in% names(cmp)))
  expect_true(all(cmp$status == "ok"))
  expect_true(all(is.finite(cmp$statistic)))
  expect_true(all(cmp$p_value >= 0 & cmp$p_value <= 1))
  expect_true(all(cmp$half_width > 0))
  expect_identical(cmp$statistic_scale,
                   c("gauge", rep("wald", 5)))
  expect_equal(attr(cmp, "level"), 0.95)
  expect_equal(unname(attr(cmp, "null.value")), c(0, 0))
  expect_output(print(cmp), "not comparable as numbers")
  expect_output(print(cmp), "Tuning actually used")
  ## fixed-b is reported before EWC.
  expect_lt(which(cmp$method == "fixedb"), which(cmp$method == "ewc"))
})

test_that("comparison rows agree with the individual test functions", {
  fit <- aersn_mean(sim_var1(200, 3, phi = 0.5, seed = 141))
  settings <- list(hac = list(kernel = "Parzen", bandwidth = "andrews"),
                   fixedb = list(b = 0.3), ewc = list(nu = 14),
                   shao = list(integration = "calendar"))
  cmp <- aersn_compare(fit, null = rep(0.1, 3), settings = settings,
                       draws = 200, seed = 2)
  for (i in seq_len(nrow(cmp))) {
    m <- cmp$method[i]
    tuning <- if (is.null(settings[[m]])) list() else settings[[m]]
    tt <- do.call(aersn_test,
                  c(list(fit, null = rep(0.1, 3), method = m), tuning,
                    list(draws = 200, seed = 2)))
    expect_equal(cmp$statistic[i], unname(tt$statistic), info = m)
    expect_equal(cmp$critical_value[i], tt$critical.value, info = m)
    expect_equal(cmp$p_value[i], tt$p.value, info = m)
    expect_identical(cmp$reject[i], tt$reject)
  }
  ## Realized tuning is reported, not just what was requested.
  expect_match(cmp$tuning[cmp$method == "hac"], "Parzen")
  expect_match(cmp$tuning[cmp$method == "hac"], "Andrews")
  expect_match(cmp$tuning[cmp$method == "fixedb"], "b_grid=0.3")
  expect_match(cmp$tuning[cmp$method == "ewc"], "nu=14")
})

test_that("a subset of methods and a linear restriction are supported", {
  fit <- aersn_mean(sim_var1(200, 3, phi = 0.4, seed = 142))
  cmp <- aersn_compare(fit, methods = c("ewc", "hull"), draws = 120, seed = 3)
  expect_identical(cmp$method, c("hull", "ewc"))   # reporting order is fixed
  A <- rbind("d12" = c(1, -1, 0))
  cmp2 <- aersn_compare(fit, null = 0, contrast = A, draws = 120, seed = 3)
  expect_equal(attr(cmp2, "q"), 1L)
  for (i in seq_len(nrow(cmp2))) {
    tt <- aersn_test(fit, null = 0, contrast = A, method = cmp2$method[i],
                     draws = 120, seed = 3)
    expect_equal(cmp2$statistic[i], unname(tt$statistic), info = cmp2$method[i])
  }
  expect_error(aersn_compare(fit, methods = c("hull", "hull")), "not repeat")
  expect_error(aersn_compare(fit, methods = "wald"), "Unknown method")
  expect_error(aersn_compare(fit, settings = list(wald = list())),
               "unknown method")
  expect_error(aersn_compare(fit, settings = list(hac = "Parzen")),
               "must be a list")
})

test_that("a method that fails is reported rather than dropped", {
  set.seed(143)
  n <- 80
  psi <- matrix(rnorm(n * 2), n, 2)
  psi <- cbind(psi, psi[, 1] - psi[, 2])         # exactly collinear
  fit <- suppressWarnings(aersn(rep(0, 3), psi))
  cmp <- aersn_compare(fit, draws = 80, seed = 4)
  expect_identical(nrow(cmp), 6L)
  expect_true(any(cmp$status != "ok"))
  failed <- cmp[cmp$status != "ok", ]
  expect_true(all(is.na(failed$statistic)))
  expect_true(all(nzchar(failed$status)))
  expect_output(print(cmp), "could not be evaluated")
})

test_that("subsetting a result table gives an object that prints", {
  fit <- aersn_mean(sim_var1(150, 2, seed = 144))
  cmp <- aersn_compare(fit, draws = 80, seed = 5)
  ct <- aersn_contrast(fit, rbind(a = c(1, -1), b = c(1, 1)), draws = 80,
                       seed = 5)
  ci <- confint(fit, draws = 80, seed = 5)
  ## Dropping columns drops the class, so the object prints as a data frame.
  sub <- cmp[, c("method", "p_value")]
  expect_false(inherits(sub, "aersn_comparison"))
  expect_s3_class(sub, "data.frame")
  expect_output(print(sub), "p_value")
  expect_null(attr(sub, "estimate"))
  sub2 <- ct[, c("estimate", "lower")]
  expect_false(inherits(sub2, "aersn_contrast"))
  expect_output(print(sub2), "estimate")
  ## Keeping every column keeps the class and the header.
  expect_s3_class(cmp[, names(cmp)], "aersn_comparison")
  expect_s3_class(cmp[1:2, ], "aersn_comparison")
  expect_output(print(cmp[1:2, ]), "Comparison of inference methods")
  expect_s3_class(ct[1, ], "aersn_contrast")
  expect_output(print(ct[1, ]), "Adjusted-range increment hull intervals")
  ## Matrix-backed intervals already drop their class on subsetting.
  expect_false(inherits(ci[, 1, drop = FALSE], "aersn_confint"))
  expect_output(print(ci[, 1, drop = FALSE]), "lower")
})
