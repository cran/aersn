## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 5)

## -----------------------------------------------------------------------------
library(aersn)
set.seed(2026)
n <- 400
e <- matrix(rnorm(2 * n), n, 2)
Y <- e
for (t in 2:n) Y[t, ] <- 0.5 * Y[t - 1, ] + e[t, ]
Y <- sweep(Y, 2, c(0.12, -0.05), `+`)

fit <- aersn_mean(Y, names = c("m1", "m2"))
fit

## -----------------------------------------------------------------------------
cmp <- aersn_compare(fit, null = c(0, 0), draws = 2000, seed = 1)
cmp

## -----------------------------------------------------------------------------
aersn_test(fit, null = c(0, 0), method = "shao", draws = 2000, seed = 1)
aersn_test(fit, null = c(0, 0), method = "hac", kernel = "Parzen",
           bandwidth = "andrews")

## -----------------------------------------------------------------------------
confint(fit, method = "ewc")
aersn_contrast(fit, rbind("m1 - m2" = c(1, -1)), method = "fixedb", b = 0.5,
               draws = 2000, seed = 1)

## ----fig.height = 3.4, fig.width = 7------------------------------------------
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(aersn_region(fit, method = "hull", draws = 2000, seed = 1),
     null = c(0, 0), main = "Increment hull")
plot(aersn_region(fit, method = "shao", draws = 2000, seed = 1),
     null = c(0, 0), main = "Quadratic self-normalization")
par(op)

## -----------------------------------------------------------------------------
grid <- expand.grid(kernel = c("Bartlett", "Parzen", "Quadratic Spectral"),
                    rule = c("short", "long", "andrews", "newey-west"),
                    stringsAsFactors = FALSE)
for (i in seq_len(nrow(grid))) {
  out <- tryCatch({
    nz <- aersn_hac_lrv(fit, kernel = grid$kernel[i],
                        bandwidth = grid$rule[i])
    sprintf("bandwidth %6.3f, lag %s", nz$tuning$bandwidth,
            nz$tuning$lag_truncation)
  }, error = function(e) "not supported")
  cat(sprintf("%-20s %-12s %s\n", grid$kernel[i], grid$rule[i], out))
}

## -----------------------------------------------------------------------------
aersn_hac_lrv(fit, kernel = "Bartlett", lag = 8)$tuning[c("bandwidth",
                                                          "lag_truncation")]

## -----------------------------------------------------------------------------
for (b in c(0.2, 0.5, 0.7, 0.9, 1)) {
  nz <- aersn_fixed_b_normalizer(fit, b = b)
  cat(sprintf("requested b = %.2f -> m = %3d, realized b = %.4f\n",
              b, nz$tuning$m, nz$tuning$b_grid))
}

## -----------------------------------------------------------------------------
max(abs(aersn_fixed_b_normalizer(fit, b = 1)$matrix -
          2 * aersn_shao_normalizer(fit)$matrix))

## -----------------------------------------------------------------------------
aersn_ewc_lrv(fit)$tuning[c("nu", "rule")]
aersn_test(fit, method = "ewc", nu = 30)$critical.value

## -----------------------------------------------------------------------------
aersn_normalizer(fit, "ldl")$notes

