test_that("the centered path has the manuscript form and zero endpoints", {
  set.seed(1)
  n <- 40; q <- 3
  psi <- matrix(rnorm(n * q), n, q)
  path <- aersn_path(psi)
  expect_s3_class(path, "aersn_path")
  expect_equal(dim(path$G), c(n + 1L, q))
  expect_equal(unname(path$G[1L, ]), rep(0, q))
  expect_equal(unname(path$G[n + 1L, ]), rep(0, q))
  C <- apply(psi, 2, cumsum)
  k <- 17
  expect_equal(unname(path$G[k + 1L, ]), (C[k, ] - (k / n) * C[n, ]) / sqrt(n))
  expect_equal(path$nodes, (0:n) / n)
  expect_equal(path$centering, "calendar")
})

test_that("a linear profile reproduces calendar centering exactly", {
  set.seed(2)
  psi <- matrix(rnorm(60), 30, 2)
  p1 <- aersn_path(psi)
  p2 <- aersn_path(psi, profile = function(r) r)
  p3 <- aersn_path(psi, profile = (0:30) / 30)
  expect_equal(p1$G, p2$G)
  expect_equal(p1$G, p3$G)
  expect_equal(p2$centering, "calendar")
})

test_that("a nonlinear profile changes the centering term as stated", {
  set.seed(3)
  n <- 25
  psi <- matrix(rnorm(n * 2), n, 2)
  tau <- ((0:n) / n)^3
  path <- aersn_path(psi, profile = tau)
  C <- apply(psi, 2, cumsum)
  k <- 10
  expect_equal(unname(path$G[k + 1L, ]),
               (C[k, ] - tau[k + 1L] * C[n, ]) / sqrt(n))
  expect_equal(path$centering, "profile")
  ## Contributions that sum to zero: the profile leaves the path unchanged.
  psi0 <- sweep(psi, 2, colMeans(psi))
  expect_equal(aersn_path(psi0, profile = tau)$G, aersn_path(psi0)$G)
})

test_that("profiles are validated", {
  expect_error(aersn_profile(c(0, 0.5, 0.9), 3), "length n \\+ 1")
  expect_error(aersn_profile(c(0.1, 0.5, 1), 2), "tau\\(0\\) = 0")
  expect_error(aersn_profile(c(0, 0.7, 0.4, 1), 3), "nondecreasing")
  expect_error(aersn_profile(c(0, NA, 1), 2), "non-finite")
  expect_error(aersn_profile("a", 2), "must be NULL")
  ok <- aersn_profile(c(0, 0.2, 0.2, 1), 3)
  expect_equal(as.numeric(ok), c(0, 0.2, 0.2, 1))
})

test_that("missing values and short series are rejected explicitly", {
  expect_error(aersn_path(c(1, NA, 3)), "missing or non-finite")
  expect_error(aersn_path(matrix(c(1, 2, Inf, 4), 2, 2)), "non-finite")
  expect_error(aersn_path(1), "At least two")
  expect_error(aersn_path("a"), "numeric")
})

test_that("parameter names are propagated", {
  psi <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("a", "b")))
  expect_equal(colnames(aersn_path(psi)$G), c("a", "b"))
  expect_equal(colnames(aersn_path(psi, names = c("x", "y"))$G), c("x", "y"))
  expect_error(aersn_path(psi, names = "x"), "length")
})

test_that("the outer-product profile has the required properties", {
  set.seed(4)
  n <- 120; p <- 3
  s <- matrix(rnorm(n * p), n, p)
  s <- sweep(s, 2, colMeans(s))
  tau <- aersn_opg_profile(s)
  expect_length(tau, n + 1L)
  expect_equal(tau[1L], 0)
  expect_equal(tau[n + 1L], 1)
  expect_true(all(diff(tau) >= 0))
  expect_lt(max(abs(tau - (0:n) / n)), 0.15)   # near linear for i.i.d. scores
  expect_true(is.finite(attr(tau, "proportionality")))
  ## trace identity: tau_k = tr(Sigma_n^{-1} Sigma_k) / p
  Sn <- crossprod(s) / n
  k <- 50
  Sk <- crossprod(s[1:k, ]) / n
  expect_equal(tau[k + 1L], sum(diag(solve(Sn, Sk))) / p)
})
