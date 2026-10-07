## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 3.5)

## -----------------------------------------------------------------------------
library(aersn)
set.seed(2026)
n <- 300
y <- 0.3 + as.numeric(arima.sim(list(ar = 0.5), n))
fit <- aersn_mean(y)
fit

## -----------------------------------------------------------------------------
plot(fit$path)

## -----------------------------------------------------------------------------
aersn_test(fit, null = 0, reference = "continuous")
aersn_test(fit, null = 0, draws = 20000, seed = 1)

## -----------------------------------------------------------------------------
aersn_test(fit, null = 0, alternative = "greater", reference = "continuous")

## -----------------------------------------------------------------------------
confint(fit, level = 0.95, draws = 20000, seed = 1)
confint(fit, level = 0.95, reference = "continuous")

## -----------------------------------------------------------------------------
ci <- confint(fit, draws = 20000, seed = 1)
aersn_gauge(fit, ci[1, "upper"])      # equals the critical value
attr(ci, "critical.value")

