## Gauge of the increment hull.
##
## Production algorithm: the dual linear program of the manuscript,
## equation (gauge-lp),
##
##   T = max_{u, l} u'z  subject to  l <= u'G_j <= l + 1,  j = 0, ..., n,
##
## which uses only the n + 1 path values.  Reference algorithms used in
## the tests: the primal program over all pairwise increments and, for
## q = 2, an exact polygon evaluation.

## Dual LP evaluation.  `G` is the (n + 1) x q centered path, `z` the
## scaled estimation error.  Returns the gauge value, the maximizing
## direction u (in the original coordinates), the offset l, and solver
## information.
gauge_lp <- function(G, z, reduce = TRUE, precondition = TRUE,
                     scale = 196L) {
  q <- ncol(G)
  z <- as.numeric(z)
  if (length(z) != q) {
    aersn_abort("`z` has length ", length(z), " but the path has ", q,
                " columns.")
  }
  if (any(!is.finite(z))) aersn_abort("`z` must be finite.")
  ranges <- unname(column_ranges(G))
  if (any(ranges <= 0)) {
    j <- which(ranges <= 0)
    aersn_abort("The centered path has zero adjusted range in coordinate(s) ",
                paste(j, collapse = ", "), "; the gauge is undefined. ",
                "This happens for constant series or contributions that ",
                "are identically zero.")
  }
  if (q == 1L) {
    value <- abs(z) / ranges
    direction <- if (z >= 0) 1 / ranges else -1 / ranges
    return(list(value = value, direction = direction,
                offset = min(direction * G), status = 0L,
                n_constraints = 2L * nrow(G), reduced = FALSE,
                preconditioned = FALSE))
  }
  ## Range scaling handles units; pivoted QR also removes near-collinearity.
  ## Gs[, pivot] = Q R, so z_Q = R^{-T} zs[pivot]. Affine equivariance
  ## preserves the gauge, and the dual direction maps back via R^{-1}.
  D <- if (precondition) 1 / ranges else rep(1, q)
  Gs <- sweep(G, 2L, D, `*`)
  zs <- z * D
  if (precondition) {
    if (nrow(Gs) < q) {
      aersn_abort("The path has too few rows for a full-dimensional hull.")
    }
    qr_path <- qr(Gs, LAPACK = TRUE)
    R <- qr.R(qr_path)
    pivot <- qr_path$pivot
    if (any(abs(diag(R)) <= max(dim(Gs)) * .Machine$double.eps *
            max(abs(R)))) {
      aersn_abort("The centered path is numerically rank deficient or too ",
                  "ill-conditioned for a reliable gauge calculation.")
    }
    Gs <- qr.Q(qr_path)
    zs <- as.numeric(backsolve(R, zs[pivot], transpose = TRUE))
  }
  reduced <- FALSE
  if (reduce && q == 2L && nrow(Gs) > 8L) {
    ## Only vertices of the convex hull of the path values can bind the
    ## range constraints (the maximum and minimum of a linear function
    ## over a point set are attained at hull vertices).
    idx <- grDevices::chull(Gs[, 1L], Gs[, 2L])
    Gs <- Gs[idx, , drop = FALSE]
    reduced <- TRUE
  }
  m <- nrow(Gs)
  ## Variables (all nonnegative for lp_solve): u = u+ - u- and the free
  ## offset l = l+ - l-.  Constraints: u'G_j - l <= 1 and -u'G_j + l <= 0.
  const <- rbind(cbind(Gs, -Gs, -1, 1), cbind(-Gs, Gs, 1, -1))
  rhs <- c(rep(1, m), rep(0, m))
  objective <- c(zs, -zs, 0, 0)
  sol <- lpSolve::lp(direction = "max", objective.in = objective,
                     const.mat = const, const.dir = rep("<=", 2L * m),
                     const.rhs = rhs, scale = scale)
  if (sol$status != 0L || !is.finite(sol$objval)) {
    msg <- switch(as.character(sol$status),
                  "2" = "infeasible",
                  "3" = "unbounded (check the path span and numerical conditioning)",
                  "5" = "numerical failure",
                  paste("lp_solve status", sol$status))
    aersn_abort("The linear program for the gauge did not solve: ", msg,
                ". Check the hull diagnostics (numerical rank, projected ",
                "ranges) of the object.")
  }
  u_scaled <- sol$solution[seq_len(q)] - sol$solution[q + seq_len(q)]
  offset <- sol$solution[2L * q + 1L] - sol$solution[2L * q + 2L]
  projection <- as.numeric(Gs %*% u_scaled) - offset
  if (any(!is.finite(projection)) ||
      min(projection) < -1e-7 || max(projection) > 1 + 1e-7) {
    aersn_abort("The gauge solver returned an infeasible direction. ",
                "Check the numerical conditioning of the centered path.")
  }
  if (precondition) {
    u_original <- numeric(q)
    u_original[pivot] <- backsolve(R, u_scaled)
    u_scaled <- u_original
  }
  direction <- u_scaled * D
  if (any(!is.finite(direction))) {
    aersn_abort("The gauge direction could not be mapped back reliably; ",
                "check the units and conditioning of the centered path.")
  }
  list(value = max(sol$objval, 0), direction = direction,
       offset = offset,
       status = 0L,
       n_constraints = 2L * m, reduced = reduced, preconditioned = precondition)
}

## Primal reference implementation (tests only): the gauge of the convex
## hull of all pairwise increments d = G_i - G_j is the value of
##   min sum(w)  subject to  sum_k w_k d_k = z,  w >= 0.
## It forms all (n + 1)^2 increments and is used on small examples.
gauge_primal <- function(G, z) {
  m <- nrow(G)
  q <- ncol(G)
  idx <- expand.grid(i = seq_len(m), j = seq_len(m))
  idx <- idx[idx$i != idx$j, , drop = FALSE]
  D <- G[idx$i, , drop = FALSE] - G[idx$j, , drop = FALSE]
  sol <- lpSolve::lp(direction = "min", objective.in = rep(1, nrow(D)),
                     const.mat = t(D), const.dir = rep("=", q),
                     const.rhs = as.numeric(z))
  if (sol$status != 0L) {
    aersn_abort("Primal reference program failed with status ", sol$status)
  }
  sol$objval
}

## Vertices of the increment hull for q = 2 (in counter-clockwise order).
## K = P - P where P is the convex hull of the path values, so the
## vertices of K are among the pairwise differences of the vertices of P.
hull_vertices_2d <- function(G) {
  if (ncol(G) != 2L) aersn_abort("hull_vertices_2d() requires q = 2.")
  idx <- grDevices::chull(G[, 1L], G[, 2L])
  V <- G[idx, , drop = FALSE]
  m <- nrow(V)
  ii <- rep(seq_len(m), times = m)
  jj <- rep(seq_len(m), each = m)
  Dm <- V[ii, , drop = FALSE] - V[jj, , drop = FALSE]
  Dm <- unique(Dm)
  idx2 <- grDevices::chull(Dm[, 1L], Dm[, 2L])
  K <- Dm[idx2, , drop = FALSE]
  ## chull() returns clockwise order; reverse to counter-clockwise.
  K <- K[rev(seq_len(nrow(K))), , drop = FALSE]
  colnames(K) <- colnames(G)
  K
}

## Exact polygon gauge for q = 2 (tests only).  With counter-clockwise
## vertices V of a convex polygon containing the origin in its interior,
## the polygon is {x: a_e'x <= b_e} over edges e, and the gauge of z is
## max_e a_e'z / b_e.
gauge_polygon_2d <- function(V, z) {
  m <- nrow(V)
  nxt <- c(seq_len(m)[-1L], 1L)
  edge <- V[nxt, , drop = FALSE] - V
  normals <- cbind(edge[, 2L], -edge[, 1L])   # outward for CCW order
  b <- rowSums(normals * V)
  if (any(b <= 0)) aersn_abort("Polygon does not contain the origin.")
  max(as.numeric(normals %*% as.numeric(z)) / b)
}
