## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 3.5)

## -----------------------------------------------------------------------------
library(aersn)
set.seed(31)
n <- 400
r <- (1:n) / n
sigma <- exp(1 - 2 * r)                    # scale falls from e^1 to e^-1
x <- rnorm(n)
beta0 <- c(0.5, 1)
X <- cbind(1, x)
y <- as.numeric(X %*% beta0) + sigma * rnorm(n)
## weighted least squares = maximum likelihood with known sigma_t
w <- 1 / sigma^2
b <- solve(crossprod(X * w, X), crossprod(X * w, y))
s <- X * as.numeric(y - X %*% b) * w      # conditional scores
J <- crossprod(X * w, X) / n               # normalized negative Hessian
tau_true <- c(0, cumsum(w)) / sum(w)

## -----------------------------------------------------------------------------
fit_cal <- aersn_mle(s, J, b, target = 2, names = "slope")
fit_opg <- aersn_mle(s, J, b, target = 2, names = "slope", profile = "opg")
plot((0:n) / n, fit_opg$nodes, type = "l", xlab = "sample fraction",
     ylab = "tau", main = "Estimated (solid) and true (dashed) profile")
lines((0:n) / n, tau_true, lty = 2)
abline(0, 1, col = "grey60", lty = 3)
fit_opg$model$proportionality      # departure from a common scalar profile

## -----------------------------------------------------------------------------
all.equal(fit_cal$path$G, fit_opg$path$G)
ref_cal <- aersn_reference(fit_cal, draws = 20000, seed = 1)
ref_opg <- aersn_reference(fit_opg, draws = 20000, seed = 1)
c(calendar = aersn_critical_value(ref_cal), profile = aersn_critical_value(ref_opg))
aersn_test(fit_opg, null = 1, reference = ref_opg)

