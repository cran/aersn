## Internal helpers: validation, deterministic directions, RNG isolation,
## and small numerical utilities.  Nothing in this file is exported.

aersn_abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

aersn_warn <- function(...) {
  warning(paste0(...), call. = FALSE)
}

## Coerce influence contributions (or any observation-by-parameter array)
## to a double matrix with observations in rows.  Missing or non-finite
## values are an error: the statistic cannot be defined without every
## contribution, and silently dropping rows would change the path.
as_contribution_matrix <- function(x, what = "psi") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    if (!is.numeric(x)) {
      aersn_abort("`", what, "` must be numeric; got ", class(x)[1L], ".")
    }
    x <- matrix(as.numeric(x), ncol = 1L)
  }
  if (length(dim(x)) != 2L) {
    aersn_abort("`", what, "` must be a numeric vector or a matrix with ",
                "observations in rows.")
  }
  if (!is.numeric(x)) {
    aersn_abort("`", what, "` must be numeric; got ", mode(x), ".")
  }
  storage.mode(x) <- "double"
  if (nrow(x) < 2L || ncol(x) < 1L) {
    aersn_abort("At least two observations and one column are required in `", what, "`.")
  }
  bad <- !is.finite(x)
  if (any(bad)) {
    rows <- which(rowSums(bad) > 0)
    aersn_abort("`", what, "` contains ", sum(bad), " missing or non-finite ",
                "value(s) in ", length(rows), " row(s) (first: ",
                paste(utils::head(rows, 5L), collapse = ", "), "). ",
                "The centered path is undefined with missing contributions; ",
                "handle missing observations before calling the package. ",
                "Rows are never dropped silently.")
  }
  x
}

check_scalar_number <- function(x, what, lower = -Inf, upper = Inf,
                                integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    aersn_abort("`", what, "` must be a single finite number.")
  }
  if (x < lower || x > upper) {
    aersn_abort("`", what, "` must lie in [", lower, ", ", upper, "]; got ",
                x, ".")
  }
  if (integer && x != round(x)) {
    aersn_abort("`", what, "` must be a whole number; got ", x, ".")
  }
  if (integer && abs(x) > .Machine$integer.max) {
    aersn_abort("`", what, "` exceeds the supported integer range.")
  }
  invisible(x)
}

check_flag <- function(x, what) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    aersn_abort("`", what, "` must be TRUE or FALSE.")
  }
  invisible(x)
}

check_probs <- function(probs) {
  if (!is.numeric(probs) || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1)) {
    aersn_abort("`probs` must lie strictly between 0 and 1.")
  }
  invisible(probs)
}

check_indices <- function(x, size, what) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
      any(x != trunc(x)) || any(x < 1 | x > size) || anyDuplicated(x)) {
    aersn_abort("`", what, "` must select distinct integer coordinates in 1..",
                size, ".")
  }
  as.integer(x)
}

check_level <- function(level) {
  check_scalar_number(level, "level", lower = 0, upper = 1)
  if (level <= 0 || level >= 1) {
    aersn_abort("`level` must lie strictly between 0 and 1; got ", level, ".")
  }
  invisible(level)
}

parameter_names <- function(names, q, estimate = NULL, psi = NULL) {
  if (!is.null(names)) {
    names <- as.character(names)
    if (length(names) != q) {
      aersn_abort("`names` has length ", length(names), " but there are ",
                  q, " parameters.")
    }
    return(names)
  }
  if (!is.null(estimate) && !is.null(names(estimate))) {
    return(names(estimate))
  }
  if (!is.null(psi) && !is.null(colnames(psi))) {
    return(colnames(psi))
  }
  if (q == 1L) "theta" else paste0("theta", seq_len(q))
}

## First `m` primes by trial division (m is at most the parameter dimension).
first_primes <- function(m) {
  primes <- integer(0)
  candidate <- 2L
  while (length(primes) < m) {
    if (all(candidate %% primes[primes <= sqrt(candidate)] != 0L)) {
      primes <- c(primes, candidate)
    }
    candidate <- candidate + 1L
  }
  primes
}

## van der Corput radical-inverse sequence in a given base.
radical_inverse <- function(index, base) {
  out <- numeric(length(index))
  denom <- 1
  k <- index
  while (any(k > 0)) {
    denom <- denom * base
    out <- out + (k %% base) / denom
    k <- k %/% base
  }
  out
}

## Deterministic, RNG-free set of unit directions in R^q (Halton points
## mapped through the normal quantile function and normalized).  Used for
## the projected-range diagnostics of the increment hull.
halton_directions <- function(m, q) {
  if (m < 1L) return(matrix(numeric(0), 0L, q))
  bases <- first_primes(q)
  idx <- seq_len(m) + 20L
  U <- vapply(bases, function(b) radical_inverse(idx, b), numeric(m))
  U <- matrix(U, nrow = m, ncol = q)
  Z <- stats::qnorm(pmin(pmax(U, 1e-10), 1 - 1e-10))
  Z / sqrt(rowSums(Z^2))
}

## Save the user's RNG state and return a function that restores it.
rng_state_saver <- function() {
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = env) else NULL
  old_kind <- RNGkind()
  function() {
    suppressWarnings(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
    invisible(NULL)
  }
}

## Deterministic seed for a simulation chunk.  Kept in double arithmetic
## and reduced modulo 2^31 - 1 so that it is a valid integer seed.
chunk_seed <- function(seed, chunk) {
  as.integer((as.double(seed) * 1000003 + 104729 * as.double(chunk)) %%
               2147483647)
}

## Short deterministic key for a numeric vector (used to identify a
## nonuniform reference grid in the cache).  Polynomial rolling hash of
## the 15-significant-digit representation modulo two large primes.
vector_key <- function(x) {
  s <- paste(formatC(x, digits = 15, format = "g"), collapse = ",")
  bytes <- utf8ToInt(s)
  h1 <- 7; h2 <- 11
  for (b in bytes) {
    h1 <- (h1 * 131 + b) %% 2147483647
    h2 <- (h2 * 137 + b) %% 2147483629
  }
  sprintf("%d-%d-%d", as.integer(h1), as.integer(h2), length(x))
}

symmetric_inverse_sqrt <- function(S, what = "matrix") {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  if (min(e$values) <= 0 ||
      min(e$values) <= 1e-12 * max(e$values)) {
    aersn_abort("`", what, "` is not numerically positive definite ",
                "(smallest eigenvalue ", signif(min(e$values), 3),
                ", largest ", signif(max(e$values), 3), ").")
  }
  e$vectors %*% (t(e$vectors) / sqrt(e$values))
}

fmt_num <- function(x, digits = 4L) {
  formatC(x, digits = digits, format = "fg", flag = "#")
}

column_ranges <- function(G) {
  apply(G, 2L, function(x) max(x) - min(x))
}
