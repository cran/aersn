## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 3.5)

## -----------------------------------------------------------------------------
library(aersn)
aersn_scalar_quantile(c(0.90, 0.95, 0.99))
curve(aersn_scalar_density(x), -4, 4, ylab = "density of M = Z / R")

## -----------------------------------------------------------------------------
sapply(c(50, 200, 1000), function(n)
  aersn_critical_value(aersn_reference(1, n = n, draws = 20000, seed = 1)))
aersn_scalar_quantile(0.95)

## -----------------------------------------------------------------------------
r1 <- aersn_reference(2, n = 100, draws = 2000, seed = 5)
r2 <- aersn_reference(2, n = 100, draws = 2000, seed = 5)
identical(r1$draws, r2$draws)
r1

## ----error = TRUE-------------------------------------------------------------
try({
set.seed(1)
fit <- aersn_mean(matrix(rnorm(300), 150, 2))
aersn_test(fit, reference = r1)
})

## -----------------------------------------------------------------------------
subset(aersn_registry, prob == 0.95)
fit200 <- aersn_mean(matrix(rnorm(400), 200, 2))
aersn_test(fit200, reference = "registry")

