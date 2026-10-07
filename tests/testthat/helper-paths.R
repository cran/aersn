## Helpers shared by the tests.

random_path <- function(n, q, seed) {
  set.seed(seed)
  aersn_path(matrix(rnorm(n * q), n, q))
}

## Brownian bridge (plus independent endpoint) on nested grids sharing
## one underlying path: returns a list of (G, z) for each grid size.
nested_bridges <- function(q, grids, seed) {
  set.seed(seed)
  m <- max(grids)
  inc <- matrix(rnorm(m * q, sd = sqrt(1 / m)), m, q)
  W <- apply(inc, 2L, cumsum)
  if (q == 1L) W <- matrix(W, ncol = 1L)
  z <- rnorm(q)
  lapply(grids, function(g) {
    step <- m %/% g
    idx <- seq.int(step, m, by = step)
    Wg <- W[idx, , drop = FALSE]
    tt <- seq_len(g) / g
    list(G = rbind(rep(0, q), Wg - tcrossprod(tt, W[m, ])), z = z)
  })
}
