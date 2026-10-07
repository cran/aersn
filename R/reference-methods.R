## Reference laws for the comparator methods.
##
## Every reference object records the statistic it belongs to and the tuning
## values that change its law, so that a statistic can never be combined with
## a reference simulated for another method, another grid, another integration
## rule, or another bandwidth fraction.

reference_statistic_label <- function(statistic) {
  switch(statistic,
         hull_gauge = "increment-hull gauge",
         componentwise_mq2 = "LDL componentwise adjusted range",
         shao_sq2 = "quadratic self-normalization",
         fixedb_bartlett = "Bartlett fixed-b",
         chisq = "chi-squared",
         ewc_F = "equal-weighted cosine (scaled F)",
         statistic)
}

## Arguments that change a simulated law, by statistic.
reference_arg_names <- function(statistic) {
  switch(statistic,
         shao_sq2 = "integration",
         fixedb_bartlett = c("b_grid", "m"),
         character())
}

check_reference_args <- function(statistic, args) {
  known <- c("hull_gauge", "componentwise_mq2", "shao_sq2", "fixedb_bartlett")
  if (!is.character(statistic) || length(statistic) != 1L ||
      !statistic %in% known) {
    aersn_abort("`statistic` must be one of ",
                paste0("\"", known, "\"", collapse = ", "), ".")
  }
  if (!is.list(args) || (length(args) &&
      (is.null(names(args)) || anyNA(names(args)) ||
       any(!nzchar(names(args))) || anyDuplicated(names(args))))) {
    aersn_abort("`args` must be a list with unique, nonempty argument names.")
  }
  need <- reference_arg_names(statistic)
  extra <- setdiff(names(args), need)
  if (length(extra)) {
    aersn_abort("Reference argument(s) ", paste(extra, collapse = ", "),
                " do not apply to statistic '", statistic, "'.")
  }
  missing <- setdiff(need, names(args))
  if (length(missing)) {
    aersn_abort("Statistic '", statistic, "' needs reference argument(s) ",
                paste(missing, collapse = ", "), ".")
  }
  if (statistic == "shao_sq2") {
    if (!is.character(args$integration) || length(args$integration) != 1L ||
        is.na(args$integration) || !args$integration %in% c("calendar", "profile")) {
      aersn_abort("`integration` must be \"calendar\" or \"profile\".")
    }
  }
  if (statistic == "fixedb_bartlett") {
    check_scalar_number(args$b_grid, "b_grid", lower = 0, upper = 1)
    check_scalar_number(args$m, "m", lower = 1, integer = TRUE)
  }
  args[need]
}

reference_args_key <- function(args) {
  if (!length(args)) return("noargs")
  paste(vapply(names(args), function(nm) {
    v <- args[[nm]]
    paste0(nm, "=", if (is.numeric(v)) formatC(v, digits = 15, format = "g")
           else as.character(v))
  }, character(1L)), collapse = "|")
}

reference_args_agree <- function(a, b) {
  a <- if (is.null(a)) list() else a
  b <- if (is.null(b)) list() else b
  if (!setequal(names(a), names(b))) return(FALSE)
  for (nm in names(a)) {
    x <- a[[nm]]; y <- b[[nm]]
    if (is.numeric(x) && is.numeric(y)) {
      if (!isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-12))) {
        return(FALSE)
      }
    } else if (!identical(x, y)) {
      return(FALSE)
    }
  }
  TRUE
}

describe_reference_args <- function(statistic, args) {
  if (is.null(args) || !length(args)) return("")
  switch(statistic,
         shao_sq2 = paste0(args$integration, "-time integration"),
         fixedb_bartlett = paste0("b = ", fmt_num(args$b_grid, 6),
                                  " (bandwidth in grid intervals ", args$m, ")"),
         paste(names(args), unlist(lapply(args, function(v) paste(v, collapse = ","))),
               sep = " = ", collapse = ", "))
}

## ---------------------------------------------------------------------------
## Matched-grid simulation for the four nonstandard laws.
## ---------------------------------------------------------------------------

simulate_matched_reference <- function(statistic, q, nodes, draws, seed, chunk,
                                       cores, args) {
  if (statistic == "hull_gauge") {
    return(simulate_gauge_reference(q, nodes, draws, seed, chunk, cores))
  }
  if (statistic == "componentwise_mq2" && q == 1L) {
    ## With one parameter the LDL transformation is the identity and the
    ## componentwise statistic is the square of the scalar adjusted-range
    ## statistic. Reusing the scalar simulation makes the two reference laws
    ## exactly consistent draw by draw.
    sim <- simulate_gauge_reference(1L, nodes, draws, seed, chunk, cores)
    return(list(values = sim$values^2, signed = NULL, failures = sim$failures))
  }
  restore <- rng_state_saver()
  on.exit(restore(), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  n <- length(nodes) - 1L
  sd_inc <- sqrt(diff(nodes))
  weights <- switch(statistic,
                    shao_sq2 = if (args$integration == "profile")
                      diff(nodes) else rep(1 / n, n),
                    NULL)
  functional <- switch(
    statistic,
    componentwise_mq2 = function(G, z) {
      r <- apply(G, 2L, function(x) max(x) - min(x))
      if (any(r <= 0)) return(NA_real_)
      sum((z / r)^2)
    },
    shao_sq2 = function(G, z) {
      V <- crossprod(G[-1L, , drop = FALSE] * sqrt(weights))
      quad_form_or_na(z, V)
    },
    fixedb_bartlett = function(G, z) {
      quad_form_or_na(z, fixed_b_matrix_from_m(G, args$m)$V)
    },
    aersn_abort("Internal: no simulator for statistic '", statistic, "'."))
  starts <- seq.int(1L, draws, by = chunk)
  sizes <- pmin(chunk, draws - starts + 1L)
  run_chunk <- function(i) {
    set.seed(chunk_seed(seed, i))
    b <- sizes[i]
    values <- numeric(b)
    failures <- 0L
    for (r in seq_len(b)) {
      inc <- matrix(stats::rnorm(n * q), nrow = n, ncol = q) * sd_inc
      W <- apply(inc, 2L, cumsum)
      if (n == 1L) W <- matrix(W, ncol = q)
      endpoint <- W[n, ]
      G <- rbind(rep(0, q), W - tcrossprod(nodes[-1L], endpoint))
      z <- stats::rnorm(q)
      val <- tryCatch(functional(G, z), error = function(e) NA_real_)
      if (!is.finite(val)) failures <- failures + 1L
      values[r] <- val
    }
    list(values = values, signed = NULL, failures = failures)
  }
  idx <- seq_along(starts)
  res <- if (cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(idx, run_chunk, mc.cores = cores, mc.preschedule = TRUE)
  } else {
    lapply(idx, run_chunk)
  }
  if (any(vapply(res, function(x) inherits(x, "try-error") || is.null(x),
                 logical(1L)))) {
    aersn_abort("A reference-simulation worker failed; no reference was returned.")
  }
  values <- unlist(lapply(res, `[[`, "values"), use.names = FALSE)
  failures <- sum(vapply(res, `[[`, integer(1L), "failures"))
  if (failures > 0L || length(values) != draws || any(!is.finite(values))) {
    aersn_abort(failures, " of ", draws, " reference draws for statistic '",
                statistic, "' failed numerically. No reference was returned ",
                "or cached; inspect the grid, the dimension and the tuning ",
                "value rather than discarding failed draws.")
  }
  list(values = values, signed = NULL, failures = failures)
}

quad_form_or_na <- function(z, V) {
  V <- (V + t(V)) / 2
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(ev)) return(NA_real_)
  drop(crossprod(z, solve(V, z)))
}

## ---------------------------------------------------------------------------
## Parametric laws: chi-squared for HAC, scaled F for EWC.
## ---------------------------------------------------------------------------

#' Parametric reference law for a Wald statistic
#'
#' Builds an [aersn_reference()] object for the two comparator methods whose
#' reference law is available in closed form: the chi-squared law used with a
#' consistent HAC covariance estimate, and the scaled F law used with the
#' equal-weighted cosine estimate.
#'
#' @details
#' For `law = "chisq"` the statistic is referred to \eqn{\chi^2_{df}}.  For
#' `law = "ewc_F"` the critical value is
#' \eqn{\nu\,df/(\nu - df + 1)} times the \eqn{F_{df,\ \nu - df + 1}} quantile,
#' equivalently \eqn{\{(\nu - df + 1)/(\nu\,df)\}\,T \sim F_{df,\nu-df+1}}.
#' Both are asymptotic reference laws for dependent observations, not exact
#' finite-sample distributions.
#'
#' @param law `"chisq"` or `"ewc_F"`.
#' @param df Degrees of freedom, that is the dimension of the parameter or of
#'   the linear restriction being tested.
#' @param nu Number of cosine terms, required for `"ewc_F"`.
#' @return An object of class `"aersn_reference"` with `type = "parametric"`.
#' @seealso [aersn_reference()], [aersn_hac_lrv()], [aersn_ewc_lrv()].
#' @examples
#' ref <- aersn_parametric_reference("chisq", df = 3)
#' aersn_critical_value(ref, 0.95)
#' aersn_parametric_reference("ewc_F", df = 2, nu = 20)
#' @export
aersn_parametric_reference <- function(law = c("chisq", "ewc_F"), df,
                                       nu = NULL) {
  law <- match.arg(law)
  check_scalar_number(df, "df", lower = 1, integer = TRUE)
  df <- as.integer(df)
  if (law == "ewc_F") {
    if (is.null(nu)) aersn_abort("`nu` is required for the EWC reference law.")
    check_scalar_number(nu, "nu", lower = 1, integer = TRUE)
    nu <- as.integer(nu)
    if (nu < df) {
      aersn_abort("The EWC reference law needs nu >= df; got nu = ", nu,
                  " and df = ", df, ".")
    }
  }
  structure(list(type = "parametric", statistic = law, law = law,
                 args = list(), q = df, df = df, nu = nu,
                 n = Inf, uniform = TRUE, nodes = NULL,
                 grid_id = "parametric", n_draws = NA_integer_,
                 seed = NA_integer_,
                 key = paste("parametric", law, "df", df,
                             if (is.null(nu)) "" else paste0("nu", nu),
                             sep = ":"),
                 version = as.character(utils::packageVersion("aersn"))),
            class = "aersn_reference")
}

parametric_quantile <- function(ref, probs) {
  if (ref$law == "chisq") return(stats::qchisq(probs, df = ref$df))
  scale <- ref$nu * ref$df / (ref$nu - ref$df + 1)
  scale * stats::qf(probs, df1 = ref$df, df2 = ref$nu - ref$df + 1)
}

parametric_tail <- function(ref, statistic) {
  if (statistic < 0) return(1)
  if (ref$law == "chisq") {
    return(stats::pchisq(statistic, df = ref$df, lower.tail = FALSE))
  }
  scale <- ref$nu * ref$df / (ref$nu - ref$df + 1)
  stats::pf(statistic / scale, df1 = ref$df, df2 = ref$nu - ref$df + 1,
            lower.tail = FALSE)
}

parametric_summary <- function(ref) {
  if (ref$law == "chisq") {
    return(sprintf("chi-squared law with %d degrees of freedom", ref$df))
  }
  sprintf(paste0("scaled F law: %d nu/(nu - %d + 1) times F(%d, %d) with ",
                 "nu = %d cosine terms"),
          ref$df, ref$df, ref$df, ref$nu - ref$df + 1, ref$nu)
}

## ---------------------------------------------------------------------------
## Resolving the reference for a normalizer.
## ---------------------------------------------------------------------------

## Statistic identifier and law-changing arguments implied by a normalizer.
normalizer_reference_spec <- function(nz, dim) {
  switch(nz$method,
         hull = list(kind = "simulated", statistic = "hull_gauge",
                     args = list()),
         ldl = list(kind = "simulated", statistic = "componentwise_mq2",
                    args = list()),
         shao = list(kind = "simulated", statistic = "shao_sq2",
                     args = list(integration = nz$tuning$integration)),
         fixedb = list(kind = "simulated", statistic = "fixedb_bartlett",
                       args = list(b_grid = nz$reference_args$b_grid,
                                   m = nz$reference_args$m)),
         hac = list(kind = "parametric", law = "chisq", df = dim),
         ewc = list(kind = "parametric", law = "ewc_F", df = dim,
                    nu = nz$tuning$nu))
}

## Build or validate the reference used with a given normalizer.
reference_for_normalizer <- function(nz, dim, reference, n, nodes, uniform,
                                     draws, seed, cores = 1L,
                                     verbose = FALSE) {
  spec <- normalizer_reference_spec(nz, dim)
  if (spec$kind == "parametric") {
    if (is.null(reference)) {
      return(aersn_parametric_reference(spec$law, df = spec$df, nu = spec$nu))
    }
    if (identical(reference, "continuous") || identical(reference, "registry")) {
      aersn_abort("Method '", nz$method, "' uses a closed-form ",
                  reference_statistic_label(spec$law), " reference law; ",
                  "reference = \"", reference, "\" does not apply to it.")
    }
    check_reference(reference)
    if (!identical(reference$type, "parametric") ||
        !identical(reference$law, spec$law) ||
        !identical(as.integer(reference$df), as.integer(spec$df)) ||
        (spec$law == "ewc_F" &&
         !identical(as.integer(reference$nu), as.integer(spec$nu)))) {
      aersn_abort("The supplied reference does not match the ",
                  reference_statistic_label(spec$law), " law required by ",
                  "method '", nz$method, "' at dimension ", spec$df,
                  if (spec$law == "ewc_F") paste0(" with nu = ", spec$nu) else "",
                  ".")
    }
    return(reference)
  }
  if (is.null(reference)) {
    ref <- aersn_reference(dim, n = n, draws = draws, seed = seed,
                           nodes = if (uniform) NULL else nodes,
                           cores = cores, verbose = verbose,
                           statistic = spec$statistic, args = spec$args)
  } else if (identical(reference, "continuous")) {
    if (spec$statistic != "hull_gauge") {
      aersn_abort("A closed-form continuous-path law is available only for ",
                  "the increment-hull gauge; method '", nz$method,
                  "' needs a matched-grid reference.")
    }
    ref <- aersn_reference(dim, grid = "continuous")
  } else if (identical(reference, "registry")) {
    if (spec$statistic != "hull_gauge") {
      aersn_abort("The bundled registry tabulates the increment-hull gauge ",
                  "only; method '", nz$method, "' needs a simulated ",
                  "reference.")
    }
    ref <- registry_reference(dim, n)
  } else {
    ref <- reference
  }
  check_reference_match(ref, q = dim, n = n, nodes = nodes, uniform = uniform,
                        statistic = spec$statistic, args = spec$args)
  ref
}
