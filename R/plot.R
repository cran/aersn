## Plot methods (base graphics).

## Let graphical arguments supplied in ... override the default labels.
plot_region_frame <- function(x, y = NULL, defaults, dots) {
  args <- c(list(x = x), if (!is.null(y)) list(y = y),
            utils::modifyList(defaults, dots))
  do.call(graphics::plot, args)
}

#' @export
plot.aersn_path <- function(x, ..., main = "Centered influence path",
                            xlab = "sample fraction", ylab = "path value") {
  r <- seq.int(0L, x$n) / x$n
  graphics::matplot(r, x$G, type = "l", lty = 1, main = main, xlab = xlab,
                    ylab = ylab, ...)
  graphics::abline(h = 0, col = "grey60", lty = 3)
  if (x$q > 1L) {
    graphics::legend("topleft", legend = x$names, col = seq_len(x$q), lty = 1,
                     bty = "n", cex = 0.8)
  }
  invisible(x)
}

#' Plot a confidence region
#'
#' For `q = 1` the interval is drawn; for `q = 2` the region polygon; for
#' `q >= 3` a matrix of two-dimensional views, either projections of the
#' joint region (`type = "projection"`, simultaneous coverage) or slices
#' through the estimate (`type = "slice"`, cross-sections).  The estimate
#' is marked, and optionally a null value.
#'
#' @param x An [aersn_region()] object.
#' @param type For `q >= 3`, `"projection"` or `"slice"`; see
#'   [aersn_projection()].
#' @param null Optional parameter vector to mark (for example a null
#'   value); it is drawn in red when outside the region.
#' @param coords For `q >= 3`, the coordinates to display (default all).
#' @param col Fill color of the region.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return The region, invisibly.
#' @export
plot.aersn_region <- function(x, type = c("projection", "slice"), null = NULL,
                              coords = NULL, col = "steelblue", ...) {
  type <- match.arg(type)
  if (!is.null(null)) {
    null <- as.numeric(null)
    if (length(null) != x$q) aersn_abort("`null` must have length q = ", x$q)
    inside <- isTRUE(aersn_contains(x, null))
  }
  fill <- grDevices::adjustcolor(col, alpha.f = 0.3)
  if (x$q == 1L) {
    ci <- confint(x)
    plot_region_frame(c(ci[1L, "lower"], ci[1L, "upper"]), c(0, 0), defaults = list(type = "n",
                   yaxt = "n", ylab = "", xlab = x$names[1L],
                   main = sprintf("%.0f%% %s confidence interval",
                                  100 * x$level, x$normalizer$label)), dots = list(...))
    graphics::segments(ci[1L, "lower"], 0, ci[1L, "upper"], 0, lwd = 4,
                       col = col)
    graphics::points(x$center, 0, pch = 19)
    if (!is.null(null)) graphics::points(null, 0, pch = 4, lwd = 2,
                                         col = if (inside) "black" else "red")
    return(invisible(x))
  }
  if (x$q == 2L) {
    V <- if (x$family == "geometric") aersn_vertices(x) else
      aersn_projection(x, c(1L, 2L))
    pts <- rbind(V, x$center, if (!is.null(null)) null)
    plot_region_frame(pts, defaults = list(type = "n", xlab = x$names[1L], ylab = x$names[2L],
                   main = sprintf("%.0f%% %s confidence region",
                                  100 * x$level, x$normalizer$label)), dots = list(...))
    graphics::polygon(V, col = fill, border = col)
    graphics::points(x$center[1L], x$center[2L], pch = 19)
    if (!is.null(null)) graphics::points(null[1L], null[2L], pch = 4, lwd = 2,
                                         col = if (inside) "black" else "red")
    return(invisible(x))
  }
  if (is.null(coords)) coords <- seq_len(x$q)
  coords <- check_indices(coords, x$q, "coords")
  k <- length(coords)
  if (k < 2L) aersn_abort("`coords` must contain at least two coordinates.")
  old <- graphics::par(mfrow = c(k - 1L, k - 1L), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (i in seq_len(k - 1L)) {
    for (j in 2:k) {
      if (j <= i) {
        graphics::plot.new()
        next
      }
      cc <- coords[c(j, i)]
      V <- if (type == "projection") aersn_projection(x, cc) else
        aersn_slice(x, cc)
      pts <- rbind(V, x$center[cc], if (!is.null(null)) null[cc])
      plot_region_frame(pts, defaults = list(type = "n", xlab = x$names[cc[1L]],
                     ylab = x$names[cc[2L]],
                     main = if (type == "projection") "projection" else "slice"),
                     dots = list(...))
      graphics::polygon(V, col = fill, border = col)
      graphics::points(x$center[cc[1L]], x$center[cc[2L]], pch = 19)
      if (!is.null(null)) {
        graphics::points(null[cc[1L]], null[cc[2L]], pch = 4, lwd = 2,
                         col = if (inside) "black" else "red")
      }
    }
  }
  invisible(x)
}

#' @export
plot.aersn <- function(x, level = 0.95, method = "hull", ..., nu = NULL,
                       reference = NULL, draws = NULL, seed = 1L) {
  dots <- list(...)
  if (!is.null(nu)) dots$nu <- nu
  tuning <- dots[intersect(names(dots), normalizer_tuning_names())]
  graphical <- dots[setdiff(names(dots), names(tuning))]
  reg <- do.call(aersn_region,
                 c(list(x, level = level, method = method), tuning,
                   list(reference = reference, draws = draws, seed = seed)))
  do.call(plot, c(list(reg), graphical))
}

## Tuning argument names across all methods, used to separate method tuning
## from graphical parameters in plot.aersn().
normalizer_tuning_names <- function() {
  c("order", "integration", "kernel", "bandwidth", "lag", "prewhite",
    "adjust", "center", "b", "nu")
}
