## Reference distributions for the increment-hull gauge.
##
## Three kinds of reference objects share the class "aersn_reference":
##   * "analytic":    q = 1, continuous path, closed-form law of |M|;
##   * "monte_carlo": matched-grid Monte Carlo draws of the gauge of the
##                    Brownian bridge on the sample grid (any q);
##   * "table":       tabulated quantiles from a verified registry.
## Every object carries the settings that determine the law, and the
## inference functions refuse to combine a statistic with a reference
## simulated for a different dimension or grid.

.aersn_cache <- new.env(parent = emptyenv())

## Fixed chunk sizes make the simulated draws reproducible for a given
## seed regardless of the number of cores used.
CHUNK_SIZE_SCALAR <- 2000L
CHUNK_SIZE_VECTOR <- 200L

aersn_default_draws <- function(q) {
  if (q == 1L) 50000L else 10000L
}

#' Reference distribution of a self-normalized statistic
#'
#' Returns the reference law of the statistic
#' \eqn{T_q = \gamma_{K(B_q)}(Z_q)}, where \eqn{B_q} is a standard
#' \eqn{q}-dimensional Brownian bridge and \eqn{Z_q} is an independent
#' standard normal vector (manuscript equation (brownian-reference)). With
#' `grid = "matched"` the bridge is evaluated on the same grid as the
#' sample path (`n` intervals, or the supplied profile nodes), which gives
#' the matched-grid reference statistic \eqn{T_q^{(n)}} whose quantiles the
#' manuscript uses as critical values.  With `grid = "continuous"` and
#' `q = 1` the closed-form law of [aersn_scalar_cdf()] is used.
#'
#' @details
#' Each Monte Carlo draw simulates Brownian increments with variances equal
#' to the node spacings, forms the bridge \eqn{W_q(\tau_j) - \tau_j W_q(1)}
#' at the nodes, and evaluates the gauge of the piecewise-linear bridge at
#' an independent standard normal vector with the same linear program used for the
#' sample statistic ([aersn_gauge()]).  For nested grids the finite-grid
#' statistic decreases to the continuous-path statistic (Supplement,
#' Lemma "Observed grids and continuous paths"), so matched-grid quantiles
#' exceed continuous-path quantiles.  Corollary "Validity of matched-grid
#' critical values" in the Supplement covers uniform and nonuniform
#' deterministic grids and, when the grid is the estimated profile,
#' conditional quantiles given the sample.
#'
#' Draws are generated in fixed-size chunks with chunk-specific seeds
#' derived from `seed`, using the Mersenne-Twister generator with inversion
#' for normal variates; the user's random-number state is saved and
#' restored.  Results are therefore reproducible for a given `seed`,
#' `draws`, `q` and grid, independently of `cores`.  Monte Carlo standard
#' errors of quantiles are estimated from `batches` interleaved batches of
#' draws.  Objects are cached in the session by a key that contains every
#' setting that changes the law (statistic, `q`, grid nodes, `draws`,
#' `seed`, chunk size, and the linear-program settings), together with
#' `batches`. The grid must contain at least `q + 1` positive increments.
#' If any draw fails numerically, simulation stops without returning or
#' caching a reference; failed draws are never discarded to estimate quantiles.
#' Matching the grid approximates the reference law and does not make
#' inference finite-sample exact for general dependent observations.
#'
#' **Computational cost.** Scalar matched-grid simulation stores several
#' matrices with `n` columns and up to 2000 rows per worker. Peak memory
#' therefore grows with `n` and `min(draws, 2000)`; decreasing `draws`
#' while keeping it above 2000 mainly reduces runtime, not peak memory.
#' For orientation, a run with `n = 20000`, 50000 draws and one worker
#' used about 3 GB of memory and one minute on an arm64 Mac; requirements
#' vary with the platform. Fewer than 2000 draws reduces peak memory but
#' increases Monte Carlo uncertainty. For `q = 1`, `grid = "continuous"`
#' avoids simulation when a continuous-path approximation is appropriate;
#' it is a different reference law, not the matched-grid calculation.
#'
#' @param x Either the parameter dimension `q`, or an [aersn()] object
#'   (in which case `q`, `n` and the profile nodes are taken from it).
#' @param n Number of grid intervals (the sample size).  Required for
#'   `grid = "matched"` unless `x` is a fit or `nodes` is supplied.
#' @param draws Number of Monte Carlo draws; defaults to 50000 for `q = 1`
#'   and 10000 otherwise.
#' @param seed Integer seed for the reference simulation.
#' @param grid `"matched"` (bridge on the sample grid) or `"continuous"`
#'   (closed-form law; available for `q = 1` only).
#' @param nodes Optional profile nodes (length `n + 1`, see
#'   [aersn_profile()]) defining a nonuniform reference grid.
#' @param batches Number of interleaved batches for Monte Carlo standard
#'   errors.
#' @param cores Number of worker processes, either 1 or 2 (forked
#'   workers on Unix-alikes; ignored on Windows).  Does not affect the
#'   draws.
#' @param cache Use and update the session cache.
#' @param verbose Print progress messages.
#' @param statistic The statistic whose law is simulated:
#'   `"hull_gauge"` (the default, the increment-hull gauge),
#'   `"componentwise_mq2"` (the LDL componentwise statistic),
#'   `"shao_sq2"` (quadratic self-normalization) or
#'   `"fixedb_bartlett"` (Bartlett fixed-b).  The chi-squared law used with
#'   HAC estimation and the scaled F law used with the equal-weighted cosine
#'   method are closed form and are built by
#'   [aersn_parametric_reference()] instead.
#' @param args A list of the tuning values that change the law of
#'   `statistic`: `integration` (`"calendar"` or `"profile"`) for
#'   `"shao_sq2"`, and `b_grid` together with the integer bandwidth `m`
#'   for `"fixedb_bartlett"`.  Empty for the other statistics.  These values
#'   enter the cache key and are checked against the statistic they are used
#'   with.
#'
#' @return An object of class `"aersn_reference"`.  For Monte Carlo
#'   references the element `draws` holds the simulated statistics in
#'   generation order (and `signed` the signed ratios when `q = 1`).
#'   Use [quantile()], [aersn_critical_value()], [aersn_pvalue()] and
#'   [aersn_mcse()] to query it.
#' @seealso [aersn_scalar_cdf()], [aersn_registry].
#' @examples
#' ref1 <- aersn_reference(1, grid = "continuous")
#' quantile(ref1, c(0.90, 0.95, 0.99))
#' ref2 <- aersn_reference(2, n = 50, draws = 400, seed = 1)
#' ref2
#' aersn_pvalue(ref2, 2.4)
#' @export
aersn_reference <- function(x, n = NULL, draws = NULL, seed = 1L,
                            grid = c("matched", "continuous"),
                            nodes = NULL, batches = 20L, cores = 1L,
                            cache = TRUE, verbose = FALSE,
                            statistic = "hull_gauge", args = list()) {
  grid <- match.arg(grid)
  check_flag(cache, "cache")
  check_flag(verbose, "verbose")
  args <- check_reference_args(statistic, args)
  if (inherits(x, "aersn")) {
    q <- x$q
    n <- x$n
    if (x$centering == "profile") nodes <- x$nodes
  } else {
    check_scalar_number(x, "q", lower = 1, integer = TRUE)
    q <- as.integer(x)
  }
  if (grid == "continuous") {
    if (statistic != "hull_gauge") {
      aersn_abort("A closed-form continuous-path law is available only for ",
                  "the increment-hull gauge. For statistic '", statistic,
                  "' use grid = \"matched\".")
    }
    if (q != 1L) {
      aersn_abort("A closed-form continuous-path law is available only for ",
                  "q = 1. For q = ", q, " use grid = \"matched\" with the ",
                  "sample size n (or a fine grid as an approximation).")
    }
    return(analytic_reference())
  }
  if (!is.null(nodes)) {
    if (!is.numeric(nodes)) aersn_abort("`nodes` must be numeric.")
    if (!is.null(n)) {
      check_scalar_number(n, "n", lower = 2, integer = TRUE)
      if (n != length(nodes) - 1L) {
        aersn_abort("`n` must equal length(nodes) - 1.")
      }
    }
    nodes <- as.numeric(nodes)
    n <- length(nodes) - 1L
    nodes <- aersn_profile(nodes, n)
  } else {
    if (is.null(n)) {
      aersn_abort("`n` (the number of grid intervals) is required for a ",
                  "matched-grid reference.")
    }
    check_scalar_number(n, "n", lower = 2, integer = TRUE)
    n <- as.integer(n)
    nodes <- aersn_profile(NULL, n)
  }
  if (statistic == "fixedb_bartlett") {
    if (args$m > n || abs(args$b_grid - args$m / n) > 1e-12) {
      aersn_abort("Fixed-b reference arguments must satisfy 1 <= m <= n ",
                  "and b_grid = m/n.")
    }
    if (!identical(attr(nodes, "centering"), "calendar")) {
      aersn_abort("A fixed-b reference requires a uniform calendar-time grid.")
    }
  }
  if (n < q + 1L) {
    aersn_abort("A matched-grid reference needs n >= q + 1 intervals.")
  }
  if (sum(diff(nodes) > 0) < q + 1L) {
    aersn_abort("A matched-grid reference needs at least q + 1 positive ",
                "grid increments; this profile gives a degenerate bridge.")
  }
  if (is.null(draws)) draws <- aersn_default_draws(q)
  check_scalar_number(draws, "draws", lower = 1, integer = TRUE)
  check_scalar_number(seed, "seed", integer = TRUE)
  check_scalar_number(batches, "batches", lower = 1, integer = TRUE)
  check_scalar_number(cores, "cores", lower = 1, upper = 2, integer = TRUE)
  draws <- as.integer(draws)
  seed <- as.integer(seed)
  uniform <- identical(attr(nodes, "centering"), "calendar")
  grid_id <- if (uniform) paste0("n", n) else
    paste0("n", n, "-nodes", vector_key(as.numeric(nodes)))
  chunk <- if (q == 1L &&
               statistic %in% c("hull_gauge", "componentwise_mq2"))
    CHUNK_SIZE_SCALAR else CHUNK_SIZE_VECTOR
  key <- paste(statistic, "q", q, grid_id, "draws", draws, "seed", seed,
               "chunk", chunk, "lp", "lpSolve-qr-scale196", "batches", batches,
               reference_args_key(args), sep = ":")
  if (cache && exists(key, envir = .aersn_cache, inherits = FALSE)) {
    cached <- get(key, envir = .aersn_cache)
    if (uniform || identical(cached$nodes, as.numeric(nodes))) return(cached)
  }
  if (verbose) {
    message(sprintf(paste0("Simulating matched-grid reference: q = %d, ",
                           "n = %d, draws = %d, seed = %d"),
                    q, n, draws, seed))
  }
  t0 <- proc.time()[["elapsed"]]
  sim <- simulate_matched_reference(statistic, q, as.numeric(nodes), draws,
                                    seed, chunk, cores, args)
  elapsed <- proc.time()[["elapsed"]] - t0
  ref <- structure(list(type = "monte_carlo", statistic = statistic,
                        args = args,
                        q = q, n = n, uniform = uniform,
                        nodes = if (uniform) NULL else as.numeric(nodes),
                        grid_id = grid_id, draws = sim$values,
                        signed = sim$signed, n_draws = draws, seed = seed,
                        chunk_size = chunk, rng = "Mersenne-Twister/Inversion",
                        batches = as.integer(batches),
                        lp_failures = sim$failures, elapsed = elapsed,
                        key = key, created = Sys.time(),
                        version = as.character(utils::packageVersion("aersn"))),
                   class = "aersn_reference")
  if (cache) assign(key, ref, envir = .aersn_cache)
  ref
}

analytic_reference <- function() {
  structure(list(type = "analytic", statistic = "hull_gauge", q = 1L,
                 n = Inf, uniform = TRUE, nodes = NULL,
                 grid_id = "continuous", n_draws = NA_integer_,
                 seed = NA_integer_, key = "hull_gauge:q:1:continuous",
                 version = as.character(utils::packageVersion("aersn"))),
            class = "aersn_reference")
}

## Simulate `draws` values of the matched-grid gauge statistic.
simulate_gauge_reference <- function(q, nodes, draws, seed, chunk, cores) {
  restore <- rng_state_saver()
  on.exit(restore(), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  n <- length(nodes) - 1L
  sd_inc <- sqrt(diff(nodes))
  starts <- seq.int(1L, draws, by = chunk)
  sizes <- pmin(chunk, draws - starts + 1L)
  run_chunk <- function(i) {
    set.seed(chunk_seed(seed, i))
    b <- sizes[i]
    if (q == 1L) {
      inc <- matrix(stats::rnorm(b * n), nrow = b, ncol = n)
      inc <- sweep(inc, 2L, sd_inc, `*`)
      W <- t(apply(inc, 1L, cumsum))
      if (n == 1L) W <- matrix(W, ncol = 1L)
      endpoint <- W[, n]
      bridge <- W - tcrossprod(endpoint, nodes[-1L])
      z <- stats::rnorm(b)
      rng_max <- pmax(0, apply(bridge, 1L, max))
      rng_min <- pmin(0, apply(bridge, 1L, min))
      signed <- z / (rng_max - rng_min)
      list(values = abs(signed), signed = signed,
           failures = as.integer(sum(!is.finite(signed))))
    } else {
      values <- numeric(b)
      failures <- 0L
      for (r in seq_len(b)) {
        inc <- matrix(stats::rnorm(n * q), nrow = n, ncol = q) * sd_inc
        W <- apply(inc, 2L, cumsum)
        if (n == 1L) W <- matrix(W, ncol = q)
        endpoint <- W[n, ]
        G <- rbind(rep(0, q), W - tcrossprod(nodes[-1L], endpoint))
        z <- stats::rnorm(q)
        val <- tryCatch(gauge_lp(G, z)$value, error = function(e) NA_real_)
        if (!is.finite(val)) failures <- failures + 1L
        values[r] <- val
      }
      list(values = values, signed = NULL, failures = failures)
    }
  }
  idx <- seq_along(starts)
  res <- if (cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(idx, run_chunk, mc.cores = cores,
                       mc.preschedule = TRUE)
  } else {
    lapply(idx, run_chunk)
  }
  if (any(vapply(res, function(x) inherits(x, "try-error") || is.null(x),
                 logical(1L)))) {
    aersn_abort("A reference-simulation worker failed; no reference was returned.")
  }
  values <- unlist(lapply(res, `[[`, "values"), use.names = FALSE)
  signed <- if (q == 1L) unlist(lapply(res, `[[`, "signed"),
                                use.names = FALSE) else NULL
  failures <- sum(vapply(res, `[[`, integer(1L), "failures"))
  if (failures > 0L || length(values) != draws || any(!is.finite(values))) {
    aersn_abort(failures, " of ", draws, " reference draws failed numerically. ",
                "No reference was returned or cached; inspect the grid and ",
                "conditioning rather than discarding failed draws.")
  }
  list(values = values, signed = signed, failures = failures)
}

#' Build a reference object from tabulated quantiles
#'
#' Wraps tabulated critical values (for example a row set of
#' [aersn_registry]) as an `"aersn_reference"` of type `"table"`.  Such a
#' reference supports [aersn_critical_value()] at the tabulated
#' probabilities only, and [aersn_pvalue()] returns the bracket of
#' tabulated tail probabilities that contains the p-value.
#'
#' @param q Parameter dimension.
#' @param n Number of grid intervals.
#' @param probs Tabulated probabilities (upper-tail quantile levels).
#' @param values Tabulated quantiles.
#' @param mcse Optional Monte Carlo standard errors of the quantiles.
#' @param source Character description of the provenance.
#' @return An `"aersn_reference"` object of type `"table"`.
#' @examples
#' reg <- subset(aersn_registry, q == 2 & n == 200)
#' ref <- aersn_table_reference(2, 200, reg$prob, reg$cv, reg$mcse,
#'                              source = "aersn_registry")
#' aersn_critical_value(ref, 0.95)
#' @export
aersn_table_reference <- function(q, n, probs, values, mcse = NULL,
                                  source = "user table") {
  check_scalar_number(q, "q", lower = 1, integer = TRUE)
  check_scalar_number(n, "n", lower = 2, integer = TRUE)
  if (n < q + 1L) aersn_abort("A table reference needs n >= q + 1 intervals.")
  check_probs(probs)
  probs <- as.numeric(probs)
  values <- as.numeric(values)
  if (length(probs) != length(values) || length(probs) < 1L) {
    aersn_abort("`probs` and `values` must have the same positive length.")
  }
  if (anyDuplicated(round(probs, 10))) {
    aersn_abort("`probs` must be distinct to ten decimal places.")
  }
  if (any(!is.finite(values)) || any(values <= 0)) {
    aersn_abort("`values` must contain finite positive quantiles.")
  }
  if (!is.null(mcse) && (!is.numeric(mcse) || length(mcse) != length(probs) ||
      any(!is.finite(mcse)) || any(mcse < 0))) {
    aersn_abort("`mcse` must contain one finite nonnegative value per quantile.")
  }
  if (!is.character(source) || length(source) != 1L || is.na(source)) {
    aersn_abort("`source` must be a single character string.")
  }
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) {
    aersn_abort("`probs` must lie strictly between 0 and 1.")
  }
  o <- order(probs)
  probs <- probs[o]
  values <- values[o]
  if (any(diff(values) < 0)) {
    aersn_abort("Tabulated quantiles must be nondecreasing in the ",
                "probability.")
  }
  if (!is.null(mcse)) mcse <- as.numeric(mcse)[o]
  structure(list(type = "table", statistic = "hull_gauge",
                 q = as.integer(q), n = as.integer(n), uniform = TRUE,
                 nodes = NULL, grid_id = paste0("n", as.integer(n)),
                 probs = probs, values = values, mcse = mcse,
                 source = source, n_draws = NA_integer_, seed = NA_integer_,
                 key = paste("hull_gauge:q", q, "table", source, sep = ":"),
                 version = as.character(utils::packageVersion("aersn"))),
            class = "aersn_reference")
}

#' Query a reference distribution
#'
#' Quantiles, critical values, p-values and Monte Carlo standard errors of
#' an [aersn_reference()] object.
#'
#' @param x,ref An `"aersn_reference"` object.
#' @param probs Probabilities (quantile levels).
#' @param level Confidence level; the critical value is the `level`
#'   quantile of the reference law.
#' @param statistic Observed statistic value(s).
#' @param alternative For `q = 1` Monte Carlo and analytic references,
#'   `"two.sided"` (default, based on \eqn{|M|}), `"greater"` or `"less"`
#'   (based on the symmetric signed-ratio law). Joint statistics are two-sided only.
#' @param type Quantile type passed to [stats::quantile()]; the default
#'   (type 8) matches the convention of the manuscript's replication code.
#' @param ... Unused.
#'
#' @return `quantile()` and `aersn_critical_value()` return numeric
#'   vectors.  `aersn_mcse()` returns the Monte Carlo standard error of each
#'   quantile estimated from interleaved batches (`NA` for analytic
#'   references).  `aersn_pvalue()` returns a list with elements `p.value`,
#'   `mcse`, `resolution` (one over the number of draws), and for table
#'   references the bracket `lower`/`upper`. Monte Carlo results also contain
#'   `conf.int`, a 95 percent binomial interval for the reference tail probability,
#'   and `conf.level`. One-sided scalar tails use symmetry: half the absolute-tail
#'   probability on one side and its complement on the other; their resolution
#'   is half that of the two-sided tail. A zero estimated tail probability is
#'   not an upper bound on the true probability. Quantile standard errors are
#'   `NA` when fewer than two batches or two draws per batch are available.
#'   Type-8 quantile interpolation and empirical tail counting can give decisions
#'   differing within one Monte Carlo resolution unit at a boundary.
#' @examples
#' ref <- aersn_reference(1, n = 100, draws = 4000, seed = 2)
#' quantile(ref, 0.95); aersn_mcse(ref, 0.95)
#' aersn_pvalue(ref, 1.9)
#' @name reference-queries
NULL

#' @rdname reference-queries
#' @export
quantile.aersn_reference <- function(x, probs = c(0.9, 0.95, 0.99),
                                     type = 8, ...) {
  check_reference(x)
  check_probs(probs)
  check_scalar_number(type, "type", lower = 1, upper = 9, integer = TRUE)
  probs <- as.numeric(probs)
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) {
    aersn_abort("`probs` must lie strictly between 0 and 1.")
  }
  out <- switch(x$type,
                analytic = aersn_scalar_quantile(probs),
                parametric = parametric_quantile(x, probs),
                monte_carlo = as.numeric(stats::quantile(x$draws, probs,
                                                         type = type,
                                                         names = FALSE)),
                table = {
                  hit <- match(round(probs, 10), round(x$probs, 10))
                  if (anyNA(hit)) {
                    aersn_abort("Tabulated reference has quantiles only at ",
                                paste(x$probs, collapse = ", "),
                                "; requested ", paste(probs, collapse = ", "),
                                ". Simulate a matched-grid reference for ",
                                "other levels.")
                  }
                  x$values[hit]
                })
  stats::setNames(out, paste0(100 * probs, "%"))
}

#' @rdname reference-queries
#' @export
aersn_critical_value <- function(ref, level = 0.95) {
  check_reference(ref)
  check_level(level)
  unname(quantile(ref, level))
}

#' @rdname reference-queries
#' @export
aersn_mcse <- function(ref, probs = c(0.9, 0.95, 0.99), type = 8) {
  check_reference(ref)
  check_probs(probs)
  check_scalar_number(type, "type", lower = 1, upper = 9, integer = TRUE)
  probs <- as.numeric(probs)
  if (ref$type %in% c("analytic", "parametric")) {
    return(stats::setNames(rep(NA_real_, length(probs)),
                           paste0(100 * probs, "%")))
  }
  if (ref$type == "table") {
    out <- if (is.null(ref$mcse)) rep(NA_real_, length(probs)) else
      ref$mcse[match(round(probs, 10), round(ref$probs, 10))]
    return(stats::setNames(out, paste0(100 * probs, "%")))
  }
  B <- ref$batches
  if (B < 2L || length(ref$draws) < 2 * B) {
    return(stats::setNames(rep(NA_real_, length(probs)), paste0(100 * probs, "%")))
  }
  id <- ((seq_along(ref$draws) - 1L) %% B) + 1L
  out <- vapply(probs, function(p) {
    qb <- vapply(split(ref$draws, id), stats::quantile, numeric(1L),
                 probs = p, type = type, names = FALSE)
    stats::sd(qb) / sqrt(length(qb))
  }, numeric(1L))
  stats::setNames(out, paste0(100 * probs, "%"))
}

#' @rdname reference-queries
#' @export
aersn_pvalue <- function(ref, statistic,
                         alternative = c("two.sided", "greater", "less")) {
  check_reference(ref)
  alternative <- match.arg(alternative)
  statistic <- as.numeric(statistic)
  if (length(statistic) != 1L || !is.finite(statistic)) {
    aersn_abort("`statistic` must be a single finite number.")
  }
  if (alternative != "two.sided" && identical(ref$type, "parametric")) {
    aersn_abort("One-sided alternatives are not defined for the ", ref$law,
                " reference law of a Wald statistic.")
  }
  if (alternative != "two.sided" && ref$q != 1L) {
    aersn_abort("One-sided alternatives are defined only for q = 1.")
  }
  if (ref$type == "parametric") {
    return(list(p.value = parametric_tail(ref, statistic), mcse = 0,
                resolution = 0, type = "parametric"))
  }
  if (ref$type == "analytic") {
    tail2 <- scalar_tail(abs(statistic))
    p <- switch(alternative,
                two.sided = tail2,
                greater = if (statistic >= 0) tail2 / 2 else 1 - tail2 / 2,
                less = if (statistic <= 0) tail2 / 2 else 1 - tail2 / 2)
    return(list(p.value = p, mcse = 0, resolution = 0, type = "analytic"))
  }
  if (ref$type == "table") {
    if (alternative != "two.sided") {
      aersn_abort("Tabulated references support two-sided p-value brackets ",
                  "only.")
    }
    alpha <- 1 - ref$probs                # tail probabilities, decreasing
    exceeds <- abs(statistic) > ref$values # TRUE at coarser levels first
    if (!any(exceeds)) {
      lower <- max(alpha); upper <- 1
    } else if (all(exceeds)) {
      lower <- 0; upper <- min(alpha)
    } else {
      k <- max(which(exceeds))
      upper <- alpha[k]; lower <- alpha[k + 1L]
    }
    return(list(p.value = upper, mcse = NA_real_, resolution = NA_real_,
                lower = lower, upper = upper, type = "table"))
  }
  R <- length(ref$draws)
  hits <- sum(ref$draws >= abs(statistic))
  tail2 <- hits / R
  ci <- c(if (hits == 0L) 0 else stats::qbeta(0.025, hits, R - hits + 1),
          if (hits == R) 1 else stats::qbeta(0.975, hits + 1, R - hits))
  p <- tail2
  se <- sqrt(tail2 * (1 - tail2) / R)
  resolution <- 1 / R
  if (alternative != "two.sided") {
    ## Use the symmetry of the scalar law for both p-values and critical
    ## values. Random sign imbalance must not give them different laws.
    p <- tail2 / 2
    ci <- ci / 2
    if ((alternative == "greater" && statistic < 0) ||
        (alternative == "less" && statistic > 0)) {
      p <- 1 - p
      ci <- 1 - rev(ci)
    }
    se <- se / 2
    resolution <- resolution / 2
    if (statistic == 0) { p <- 0.5; ci <- c(0.5, 0.5); se <- 0 }
  }
  list(p.value = p, mcse = se, resolution = resolution,
       conf.int = ci, conf.level = 0.95, type = "monte_carlo")
}

check_reference <- function(ref) {
  if (!inherits(ref, "aersn_reference")) {
    aersn_abort("`ref` must be an aersn_reference object.")
  }
  if (!is.character(ref$type) || length(ref$type) != 1L ||
      is.na(ref$type) ||
      !ref$type %in% c("analytic", "parametric", "monte_carlo", "table")) {
    aersn_abort("Unknown reference type; rebuild the reference with aersn_reference().")
  }
  if (ref$type == "monte_carlo" &&
      (!is.numeric(ref$draws) || !length(ref$draws) ||
       any(!is.finite(ref$draws)) || any(ref$draws < 0) ||
       !isTRUE(ref$n_draws == length(ref$draws)) ||
       !isTRUE(ref$lp_failures == 0L))) {
    aersn_abort("The Monte Carlo reference is incomplete or contains failed draws; ",
                "simulate a new reference instead of using this object.")
  }
  invisible(ref)
}

## Ensure that a reference object matches the statistic it is used for.
check_reference_match <- function(ref, q, n, nodes, uniform,
                                  statistic = "hull_gauge", args = list()) {
  check_reference(ref)
  if (!identical(ref$statistic, statistic)) {
    aersn_abort("The reference object is for statistic '", ref$statistic,
                "', but the test statistic is '", statistic,
                "'. Reference laws are not interchangeable across methods.")
  }
  if (!reference_args_agree(ref$args, args)) {
    aersn_abort("The reference object was built with ",
                describe_reference_args(ref$statistic, ref$args),
                ", but the statistic uses ",
                describe_reference_args(statistic, args),
                ". Simulate a reference for the tuning value in use.")
  }
  if (ref$q != q) {
    aersn_abort("The reference distribution was computed for q = ", ref$q,
                " but the statistic has q = ", q, ".")
  }
  if (ref$type == "analytic") return(invisible(ref))
  if (ref$n != n) {
    aersn_abort("The reference distribution was computed on a grid with n = ",
                ref$n, " intervals but the sample has n = ", n, ".")
  }
  if (!uniform) {
    if (ref$uniform || !isTRUE(all.equal(ref$nodes, as.numeric(nodes),
                                         tolerance = 1e-12))) {
      aersn_abort("The statistic uses profile centering, but the reference ",
                  "distribution was not simulated on the same profile ",
                  "nodes. Simulate the reference from the fitted object ",
                  "(aersn_reference(fit)).")
    }
  } else if (!ref$uniform) {
    aersn_abort("The reference distribution was simulated on a nonuniform ",
                "profile grid, but the statistic uses calendar-time ",
                "centering.")
  }
  invisible(ref)
}

## Resolve the reference used by inference functions.
resolve_reference <- function(reference, q, n, nodes, uniform, draws, seed,
                              cores = 1L, verbose = FALSE) {
  if (is.null(reference)) {
    ref <- aersn_reference(q, n = n, draws = draws, seed = seed,
                           nodes = if (uniform) NULL else nodes,
                           cores = cores, verbose = verbose)
  } else if (identical(reference, "continuous")) {
    ref <- aersn_reference(q, grid = "continuous")
  } else if (identical(reference, "registry")) {
    ref <- registry_reference(q, n)
  } else {
    ref <- reference
  }
  check_reference_match(ref, q, n, nodes, uniform)
  ref
}

registry_reference <- function(q, n) {
  reg <- get_registry()
  hit <- reg[reg$q == q & reg$n == n, , drop = FALSE]
  if (nrow(hit) == 0L) {
    aersn_abort("The bundled registry has no entry for q = ", q, " and n = ",
                n, ". Available (q, n) pairs: ",
                paste(unique(paste0("(", reg$q, ", ", reg$n, ")")),
                      collapse = " "), ". Use a simulated reference instead.")
  }
  aersn_table_reference(q, n, hit$prob, hit$cv, hit$mcse,
                        source = paste0("aersn_registry (", hit$source[1L],
                                        ")"))
}

get_registry <- function() {
  env <- new.env()
  utils::data("aersn_registry", package = "aersn", envir = env)
  get("aersn_registry", envir = env)
}

#' Clear the session cache of reference distributions
#'
#' @return Invisibly, the number of cached objects removed.
#' @export
aersn_clear_cache <- function() {
  keys <- ls(envir = .aersn_cache, all.names = TRUE)
  rm(list = keys, envir = .aersn_cache)
  invisible(length(keys))
}

reference_summary <- function(ref) {
  extra <- describe_reference_args(ref$statistic, ref$args)
  extra <- if (nzchar(extra)) paste0(", ", extra) else ""
  switch(ref$type,
         analytic = "closed-form continuous-path law of |M| (q = 1)",
         parametric = parametric_summary(ref),
         monte_carlo = sprintf(paste0("matched-grid Monte Carlo law for %s: ",
                                      "q = %d, %s grid with n = %d intervals, ",
                                      "%d draws, seed %d%s"),
                               reference_statistic_label(ref$statistic),
                               ref$q, if (ref$uniform) "uniform" else "profile",
                               ref$n, ref$n_draws, ref$seed, extra),
         table = sprintf("tabulated quantiles for %s: q = %d, n = %d (%s)",
                         reference_statistic_label(ref$statistic),
                         ref$q, ref$n, ref$source))
}

#' @export
print.aersn_reference <- function(x, ...) {
  cat("Reference distribution (aersn_reference): ",
      reference_statistic_label(x$statistic), "\n", sep = "")
  cat("  ", reference_summary(x), "\n", sep = "")
  probs <- if (x$type == "table") x$probs else c(0.9, 0.95, 0.99)
  qs <- quantile(x, probs)
  se <- aersn_mcse(x, probs)
  cat("  quantiles (type 8):\n")
  for (i in seq_along(probs)) {
    cat(sprintf("    %5.1f%%: %s", 100 * probs[i], fmt_num(qs[i], 5)))
    if (is.finite(se[i])) cat(sprintf("  (Monte Carlo s.e. %s)",
                                      fmt_num(se[i], 2)))
    cat("\n")
  }
  invisible(x)
}
