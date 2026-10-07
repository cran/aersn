test_that("q = 1 reduces exactly to the adjusted-range ratio", {
  path <- random_path(50, 1, seed = 5)
  G <- path$G
  for (z in c(-2, -0.3, 0, 0.7, 3)) {
    g <- aersn:::gauge_lp(G, z)
    expect_equal(g$value, abs(z) / (max(G) - min(G)))
    expect_equal(g$direction * (max(G) - min(G)), sign(z) + (z == 0))
  }
})

test_that("the dual program agrees with the primal program over all increments", {
  for (q in c(2L, 3L, 5L)) {
    for (rep in 1:4) {
      path <- random_path(if (q == 5L) 9L else 12L, q, seed = 100 * q + rep)
      z <- rnorm(q)
      dual <- aersn:::gauge_lp(path$G, z)$value
      primal <- aersn:::gauge_primal(path$G, z)
      expect_equal(dual, primal, tolerance = 1e-9)
      ## the reduced (q = 2) and unreduced programs coincide
      expect_equal(aersn:::gauge_lp(path$G, z, reduce = FALSE)$value, dual,
                   tolerance = 1e-9)
      expect_equal(aersn:::gauge_lp(path$G, z, precondition = FALSE)$value,
                   dual, tolerance = 1e-9)
    }
  }
})

test_that("for q = 2 the program agrees with the exact polygon gauge", {
  for (rep in 1:5) {
    path <- random_path(60, 2, seed = 200 + rep)
    z <- rnorm(2)
    V <- aersn:::hull_vertices_2d(path$G)
    expect_equal(aersn:::gauge_lp(path$G, z)$value,
                 aersn:::gauge_polygon_2d(V, z), tolerance = 1e-9)
    ## every vertex of K has gauge one; every path increment is inside
    expect_equal(apply(V, 1, function(v) aersn:::gauge_lp(path$G, v)$value),
                 rep(1, nrow(V)), tolerance = 1e-8)
    inc <- path$G[7, ] - path$G[30, ]
    expect_lte(aersn:::gauge_lp(path$G, inc)$value, 1 + 1e-8)
  }
})

test_that("the gauge is affine equivariant, symmetric and homogeneous", {
  for (q in c(2L, 3L, 5L, 10L)) {
    path <- random_path(80, q, seed = 300 + q)
    G <- path$G
    z <- rnorm(q)
    base <- aersn:::gauge_lp(G, z)$value
    set.seed(400 + q)
    Q <- qr.Q(qr(matrix(rnorm(q * q), q, q)))
    D <- diag(exp(seq(-1.5, 1.5, length.out = q)))
    S <- diag(q); S[1, 2] <- 0.8; if (q > 2) S[3, 1] <- -0.6
    H <- matrix(rnorm(q * q), q, q)
    for (A in list(Q, D, S, H)) {
      expect_equal(aersn:::gauge_lp(G %*% t(A), as.numeric(A %*% z))$value,
                   base, tolerance = 1e-8)
    }
    expect_equal(aersn:::gauge_lp(G, -z)$value, base, tolerance = 1e-9)
    expect_equal(aersn:::gauge_lp(G, 2.5 * z)$value, 2.5 * base,
                 tolerance = 1e-9)
    expect_equal(aersn:::gauge_lp(G, 0 * z)$value, 0)
    ## translation of the path leaves the statistic unchanged
    shift <- matrix(rnorm(q), nrow(G), q, byrow = TRUE)
    expect_equal(aersn:::gauge_lp(G + shift, z)$value, base, tolerance = 1e-8)
  }
})

test_that("support function equals the projected adjusted range and bounds the gauge", {
  path <- random_path(40, 3, seed = 7)
  hull <- aersn_hull(path)
  set.seed(8)
  U <- matrix(rnorm(30), 10, 3)
  h <- aersn_support(hull, U)
  G <- path$G
  m <- nrow(G)
  for (i in seq_len(nrow(U))) {
    proj <- as.numeric(G %*% U[i, ])
    expect_equal(h[i], max(proj) - min(proj))
    ## max over all pairwise increments
    idx <- expand.grid(a = 1:m, b = 1:m)
    inc <- proj[idx$a] - proj[idx$b]
    expect_equal(h[i], max(inc))
  }
  z <- c(0.4, -1.1, 0.9)
  det <- aersn_gauge(hull, z, details = TRUE)
  ## dual representation: gauge >= |u'z| / h(u) for every u, equality at u*
  ratios <- abs(as.numeric(U %*% z)) / h
  expect_true(all(ratios <= det$value + 1e-9))
  ustar <- det$direction
  expect_equal(abs(sum(ustar * z)) / aersn_support(hull, ustar), det$value,
               tolerance = 1e-8)
  expect_equal(aersn_support(hull, ustar), 1, tolerance = 1e-8)
  ## triangle inequality
  z2 <- c(-0.3, 0.2, 0.5)
  expect_lte(aersn_gauge(hull, z + z2),
             aersn_gauge(hull, z) + aersn_gauge(hull, z2) + 1e-9)
})

test_that("finite-grid gauges decrease along nested grids sharing a path", {
  for (q in c(1L, 2L, 3L)) {
    for (rep in 1:5) {
      br <- nested_bridges(q, c(50L, 100L, 200L, 400L), seed = 500 + 10 * q + rep)
      g <- vapply(br, function(b) aersn:::gauge_lp(b$G, b$z)$value, numeric(1))
      expect_true(all(diff(g) <= 1e-9))
    }
  }
})
