## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(aersn)
set.seed(3)
n <- 200
## A ratio of two means, with the contributions derived by the delta method.
Y <- cbind(rnorm(n, 2), rnorm(n, 1))
m <- colMeans(Y)
ratio <- m[2] / m[1]
J <- matrix(c(-m[2] / m[1]^2, 1 / m[1]), 1, 2)   # derivative of h(m) = m2 / m1
psi <- sweep(Y, 2, m) %*% t(J)                     # n x 1 contributions
fit <- aersn(ratio, psi, n = n, names = "ratio")
fit
aersn_test(fit, null = 0.5, draws = 20000, seed = 1)

## -----------------------------------------------------------------------------
fit2 <- aersn_target(aersn_mean(Y), function(b) b[2] / b[1], jacobian = J,
                     names = "ratio")
all.equal(fit$psi, fit2$psi, check.attributes = FALSE)

## -----------------------------------------------------------------------------
x <- rnorm(n) + 1
y <- 1 + 0.5 * x + rnorm(n)
Z <- cbind(1, x); b <- solve(crossprod(Z), crossprod(Z, y)); u <- as.numeric(y - Z %*% b)
correct <- ((Z * u) %*% solve(crossprod(Z) / n))[, 2]
naive <- x * u / mean(x^2)
c(correct = aersn_gauge(aersn(b[2], correct), 0.5),
  naive = aersn_gauge(aersn(b[2], naive), 0.5))

