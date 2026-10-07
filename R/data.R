#' Verified matched-grid critical values from the manuscript replication
#' package
#'
#' Tabulated matched-grid quantiles of the increment-hull gauge statistic
#' \eqn{T_q^{(n)}} taken from the replication package of Hong, Lin, Linton,
#' Newey and Sun (2026), version of 8 September 2026.  They are included as
#' regression targets for the package's own simulations and as a fast
#' tabulated reference ([aersn_test()] with `reference = "registry"`).
#'
#' @format A data frame with one row per (dimension, grid size,
#'   probability):
#' \describe{
#'   \item{q}{parameter dimension (1, 2, 3, 5, 10)}
#'   \item{n}{number of grid intervals (200, 500, 738, 1000)}
#'   \item{prob}{quantile level (0.90, 0.95, 0.975, 0.99)}
#'   \item{cv}{matched-grid quantile (type-8 sample quantile of the draws)}
#'   \item{mcse}{Monte Carlo standard error from 100 interleaved batches}
#'   \item{draws}{number of Monte Carlo draws (2,000,000 for `q = 1`,
#'     30,000 otherwise)}
#'   \item{batches}{number of batches used for `mcse`}
#'   \item{quantile_type}{quantile type (8)}
#'   \item{seed}{master seed recorded in the production run manifest}
#'   \item{source}{file of the replication package the row was taken from}
#' }
#' @details The `q = 1` rows were produced by
#'   `scalar_matched_range_quadratic_cv.R` (endpoint an independent
#'   standard normal; bridge the equispaced linear interpolant; range
#'   including the zero at the origin).  The `q >= 2` rows were produced by
#'   `multivariate_matched_cv.R` with the same dual linear program used by
#'   [aersn_gauge()].  The production run manifests record master seeds
#'   20260715 and 20260717 and the replication counts stored in the
#'   `draws` column (runs completed 15 July 2026). These tabulated simulation
#'   results are distributed under this package's MIT license.
#' @source Replication package, files
#'   `core_methods/inputs/matched_scalar_critical_values.csv` and
#'   `core_methods/inputs/multivariate_matched_cv_registry.csv`.
#' @examples
#' subset(aersn_registry, prob == 0.95)
"aersn_registry"
