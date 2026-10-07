## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 5)

## -----------------------------------------------------------------------------
library(aersn)
set.seed(21)
n <- 300
x1 <- as.numeric(arima.sim(list(ar = 0.6), n))
x2 <- as.numeric(arima.sim(list(ar = 0.3), n))
u <- as.numeric(arima.sim(list(ar = 0.5), n))
y <- 1 + 0.8 * x1 - 0.4 * x2 + u
fit <- aersn_lm(y ~ x1 + x2, target = c("x1", "x2"))
fit
ref <- aersn_reference(fit, draws = 4000, seed = 1)
aersn_test(fit, null = c(0.8, -0.4), reference = ref)
confint(fit, reference = ref)
plot(fit, reference = ref, null = c(0.8, -0.4))

## -----------------------------------------------------------------------------
set.seed(22)
z1 <- rnorm(n); z2 <- rnorm(n); v <- rnorm(n)
x <- z1 + 0.5 * z2 + v
y <- 1 + x + 0.6 * v + rnorm(n)
Z <- cbind(1, z1, z2); X <- cbind(1, x)
W <- solve(crossprod(Z) / n)
Dm <- -crossprod(Z, X) / n
gZy <- crossprod(Z, y) / n
beta <- -solve(t(Dm) %*% W %*% Dm, t(Dm) %*% W %*% gZy)   # 2SLS
moments <- Z * as.numeric(y - X %*% beta)
gmm <- aersn_gmm(moments, Dm, beta, weight = W, target = 2, names = "slope")
gmm
aersn_test(gmm, null = 1, draws = 20000, seed = 1)

