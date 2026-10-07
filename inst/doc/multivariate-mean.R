## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 5)

## -----------------------------------------------------------------------------
library(aersn)
set.seed(11)
n <- 250; q <- 3
A <- matrix(c(0.4, 0.1, 0, -0.2, 0.3, 0.1, 0, 0.1, 0.5), 3, 3, byrow = TRUE)
e <- matrix(rnorm(n * q), n, q) %*% chol(matrix(c(1, .5, .2, .5, 1, .3, .2, .3, 1), 3))
Y <- e
for (t in 2:n) Y[t, ] <- Y[t - 1, ] %*% t(A) + e[t, ]
Y <- sweep(Y, 2, c(0, 0.2, -0.1), `+`)
fit <- aersn_mean(Y, names = c("m1", "m2", "m3"))
fit

## -----------------------------------------------------------------------------
ref <- aersn_reference(fit, draws = 4000, seed = 1)   # matched grid, q = 3, n = 250
ref
aersn_test(fit, null = c(0, 0, 0), reference = ref)
reg <- aersn_region(fit, reference = ref)
reg

## -----------------------------------------------------------------------------
aersn_contains(reg, rbind(c(0, 0, 0), c(0, 0.2, -0.1)))

## -----------------------------------------------------------------------------
plot(reg, null = c(0, 0, 0))

## -----------------------------------------------------------------------------
contrasts <- rbind("m2 - m1" = c(-1, 1, 0),
                   "m3 - m2" = c(0, -1, 1),
                   "average" = c(1, 1, 1) / 3)
aersn_contrast(fit, contrasts, reference = ref)

## -----------------------------------------------------------------------------
confint(fit, reference = ref)
confint(fit, type = "marginal", draws = 20000, seed = 1)

## -----------------------------------------------------------------------------
H <- matrix(c(2, 0.5, 0, -0.3, 1, 0.2, 0, 0, 0.7), 3, 3)
fitH <- aersn_mean(Y %*% t(H))
c(original = aersn_gauge(fit, c(0, 0, 0)), transformed = aersn_gauge(fitH, c(0, 0, 0)))

