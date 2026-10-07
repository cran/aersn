test_that("the sample-mean interface is the direct construction", {
  set.seed(20)
  Y <- matrix(rnorm(90), 30, 3)
  fit <- aersn_mean(Y)
  direct <- aersn(colMeans(Y), sweep(Y, 2, colMeans(Y)))
  expect_equal(fit$path$G, direct$path$G)
  expect_equal(fit$estimate, direct$estimate)
  expect_equal(fit$model$type, "mean")
  fit1 <- aersn_mean(Y[, 1])
  expect_equal(fit1$q, 1L)
})

test_that("OLS contributions equal the closed form and the GMM route", {
  set.seed(21)
  n <- 100
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  fit <- aersn_lm(y ~ x1 + x2)
  X <- cbind(1, x1, x2)
  b <- solve(crossprod(X), crossprod(X, y))
  u <- as.numeric(y - X %*% b)
  psi <- (X * u) %*% solve(crossprod(X) / n)
  expect_equal(unname(fit$psi), unname(psi))
  expect_equal(unname(fit$estimate), as.numeric(b))
  expect_equal(fit$q, 3L)
  ## lm object input, target selection by name and index
  lmfit <- lm(y ~ x1 + x2)
  fit2 <- aersn_lm(lmfit, target = "x1")
  fit3 <- aersn_lm(y ~ x1 + x2, target = 2)
  expect_equal(fit2$psi, fit3$psi)
  expect_equal(unname(fit2$psi[, 1]), psi[, 2])
  ## GMM with identity weights on the same moments
  g <- aersn_gmm(X * u, -crossprod(X) / n, b)
  expect_equal(unname(g$psi), unname(psi))
  ## moments evaluated away from the solution violate the first-order
  ## condition and trigger a warning
  u_wrong <- as.numeric(y - X %*% (b + 0.5))
  expect_warning(aersn_gmm(X * u_wrong, -crossprod(X) / n, b + 0.5),
                 "first-order condition")
  expect_error(aersn_gmm(X * u, matrix(1, 3, 4), rep(0, 4)), "Under-identified")
  expect_error(aersn_gmm(X * u, cbind(c(1, 1, 1), c(2, 2, 2)), c(1, 1)),
               "full column rank")
  expect_error(aersn_lm(lm(y ~ x1, weights = rep(1, n))), "Weighted")
  expect_error(aersn_lm(lm(y ~ x1 + I(2 * x1))), "aliased")
  expect_error(aersn_lm(3), "formula or an lm")
})

test_that("GMM warning allows accurate numerical solutions and detects residual errors", {
  set.seed(502)
  n <- 150
  x <- rnorm(n)
  y <- rpois(n, exp(0.3 + 0.5 * x))
  X <- cbind(1, x)
  moments <- function(b) X * as.numeric(y - exp(X %*% b))
  objective <- function(b) sum(colMeans(moments(b))^2)
  b <- optim(c(0, 0), objective, method = "BFGS")$par
  D <- -crossprod(X * as.numeric(exp(X %*% b)), X) / n
  expect_no_warning(fit <- aersn_gmm(moments(b), D, b))
  expect_lt(max(abs(fit$model$foc_standardized)), fit$model$foc_tolerance)
  exact <- coef(glm(y ~ x, family = poisson,
                    control = glm.control(epsilon = 1e-14)))
  De <- -crossprod(X * as.numeric(exp(X %*% exact)), X) / n
  exact_fit <- aersn_gmm(moments(exact), De, exact)
  expect_equal(aersn_gauge(fit, c(0.3, 0.5)),
               aersn_gauge(exact_fit, c(0.3, 0.5)), tolerance = 1e-4)
  S <- diag(c(3, 0.2))
  expect_no_warning(scaled <- aersn_gmm(moments(b) %*% S, S %*% D, b,
                                        weight = solve(S %*% S)))
  expect_equal(scaled$model$foc_standardized, fit$model$foc_standardized,
               tolerance = 1e-8)
  expect_warning(aersn_gmm(moments(b + 0.05), D, b + 0.05),
                 "Check optimizer convergence")
})

test_that("GMM warning uses the sqrt(n) scale and records its numerical tolerance", {
  for (n in c(40, 400)) {
    set.seed(n)
    y <- rnorm(n)
    tolerance <- .Machine$double.eps^0.25
    b <- mean(y) + 0.5 * tolerance * sd(y) / sqrt(n)
    expect_no_warning(fit <- aersn_gmm(matrix(y - b), matrix(-1), b))
    expect_equal(abs(fit$model$foc_standardized), 0.5 * tolerance,
                 tolerance = 1e-8)
    b <- mean(y) + 2 * tolerance * sd(y) / sqrt(n)
    expect_warning(aersn_gmm(matrix(y - b), matrix(-1), b),
                   "first-order condition")
  }
})

test_that("two-stage least squares is handled as over-identified GMM", {
  set.seed(22)
  n <- 200
  z1 <- rnorm(n); z2 <- rnorm(n); v <- rnorm(n)
  x <- z1 + 0.5 * z2 + v
  y <- 1 + x + 0.6 * v + rnorm(n)
  Z <- cbind(1, z1, z2); X <- cbind(1, x)
  W <- solve(crossprod(Z) / n)
  D <- -crossprod(Z, X) / n
  b <- -solve(t(D) %*% W %*% D, t(D) %*% W %*% (crossprod(Z, y) / n))
  mom <- Z * as.numeric(y - X %*% b)
  fit <- aersn_gmm(mom, D, b, weight = W, target = 2, names = "slope")
  ## closed form: psi_t = -h M m_t with M = (D'WD)^{-1} D'W
  M <- solve(t(D) %*% W %*% D, t(D) %*% W)
  psi <- -(mom %*% t(M))[, 2]
  expect_equal(unname(fit$psi[, 1]), unname(psi))
  expect_equal(max(abs(fit$model$foc)), 0, tolerance = 1e-8)
  expect_error(aersn_gmm(mom, D, b, weight = W[1:2, 1:2]), "3 x 3")
  Wbad <- W; Wbad[1, 2] <- Wbad[1, 2] + 1
  expect_error(aersn_gmm(mom, D, b, weight = Wbad), "symmetric")
})

test_that("the score interface reproduces OLS for the Gaussian linear model", {
  set.seed(23)
  n <- 150
  x <- rnorm(n); y <- 0.5 + x + rnorm(n)
  X <- cbind(1, x)
  b <- solve(crossprod(X), crossprod(X, y))
  s <- X * as.numeric(y - X %*% b)
  J <- crossprod(X) / n
  fm <- aersn_mle(s, J, b, target = 2)
  fl <- aersn_lm(y ~ x, target = 2)
  expect_equal(unname(fm$psi), unname(fl$psi))
  expect_equal(unname(fm$estimate), unname(fl$estimate))
  ## profile: fitted scores sum to zero, so the path is unchanged
  fo <- aersn_mle(s, J, b, target = 2, profile = "opg")
  expect_equal(fo$path$G, fm$path$G)
  expect_equal(fo$centering, "profile")
  expect_true(all(diff(fo$nodes) >= 0))
  expect_equal(range(fo$nodes), c(0, 1))
  fk <- aersn_mle(s, J, b, target = 2, profile = ((0:n) / n)^2)
  expect_equal(fk$nodes, ((0:n) / n)^2)
  expect_warning(aersn_mle(s + 1, J, b), "fitted scores")
  expect_error(aersn_mle(s, matrix(0, 2, 2), b), "singular")
  expect_error(aersn_mle(s, J, 1), "length p")
})

test_that("target transformations handle nuisance parameters through the Jacobian", {
  set.seed(24)
  n <- 100
  x <- rnorm(n); y <- 2 + 0.5 * x + rnorm(n)
  X <- cbind(1, x)
  b <- solve(crossprod(X), crossprod(X, y))
  u <- as.numeric(y - X %*% b)
  h <- list(h = function(beta) beta[2] / beta[1],
            jacobian = function(beta) matrix(c(-beta[2] / beta[1]^2,
                                               1 / beta[1]), 1, 2))
  fit <- aersn_gmm(X * u, -crossprod(X) / n, b, target = h, names = "ratio")
  psi_beta <- (X * u) %*% solve(crossprod(X) / n)
  J <- h$jacobian(b)
  expect_equal(unname(fit$psi[, 1]), as.numeric(psi_beta %*% t(J)))
  expect_equal(unname(fit$estimate), b[2] / b[1])
  ## the same result through aersn_target on the full fit
  full <- aersn_lm(y ~ x)
  tf <- aersn_target(full, h$h, jacobian = h$jacobian, names = "ratio")
  expect_equal(tf$psi, fit$psi, ignore_attr = TRUE)
  expect_error(aersn_gmm(X * u, -crossprod(X) / n, b, target = "zz"),
               "Unknown parameter")
  expect_error(aersn_gmm(X * u, -crossprod(X) / n, b, target = list(h = 1)),
               "function")
})
