# Leave one core free so the session stays responsive. availableCores() knows
# about cgroup and scheduler limits, which detectCores() does not, and returns
# 2 under R CMD check.
default_cores <- function() {
  max(1L, parallelly::availableCores() - 1L)
}

#' Bootstrap confidence intervals for local average treatment effects
#'
#' Resamples the data with replacement, re-estimates `P_Z`, `Q_Z`, `P_Sigma`
#' and the effects on each resample, and takes empirical quantiles of the
#' resulting distribution.
#'
#' @param n Number of bootstrap replicates. With multiply imputed data this is
#'   the number drawn from *each* imputation.
#' @param data A data frame, or a list of data frames holding multiple
#'   imputations of the same study. With a list, `n` replicates are drawn from
#'   each imputation and all replicates are pooled before taking quantiles.
#' @param instrument,treatment,outcome Column names, as strings.
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_pairs()].
#' @param alpha Two-sided error rate; the interval runs from the `alpha / 2` to
#'   the `1 - alpha / 2` quantile.
#' @param n_cores Number of worker processes. Defaults to one fewer than the
#'   cores available to this session. This affects only how long the run takes:
#'   each replicate draws from its own random number stream, so the estimates
#'   are the same whatever `n_cores` is set to.
#' @param cap If `TRUE`, discard replicates whose estimated effect lies outside
#'   the observed range of the outcome before taking quantiles. Ratio
#'   estimators with a small denominator occasionally produce extreme values,
#'   and those replicates would otherwise dominate the quantiles.
#' @param covariates,parametric,family Passed to [estimate_p_z()] and
#'   [estimate_q_z()] to bootstrap the covariate-adjusted estimator.
#' @param seed Seed for the parallel random number streams.
#' @param instrument_levels,treatments Passed to the estimators to fix the row
#'   and column order. Default to the dimnames of the projection matrices.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{`alpha`}{the error rate used.}
#'     \item{`replicates`}{a data frame with one row per bootstrap replicate
#'       and one column per identified effect.}
#'     \item{`ci`}{a data frame of `lower` and `upper` bounds, plus the number
#'       of replicates each interval is based on.}
#'   }
#'
#' @details
#' Replicates in which an instrument value drops out of the resample fail to
#' estimate, and are recorded as `NA` rather than aborting the run. The `n_kept`
#' column of `ci` reports how many replicates each interval actually used, so
#' that a badly behaved contrast is visible rather than silent.
#'
#' Given the same `seed`, the same data and the same number of replicates, the
#' result does not depend on `n_cores` or on the machine it runs on.
#'
#' Called `BSCICalculator()` in the paper code, where the multiple-imputation
#' case was selected with `Data.complete = FALSE`.
#'
#' @examples
#' \donttest{
#' treatments <- c("a", "b")
#' instrument <- list(`1` = c("a", "b"), `2` = c("b", "a"))
#' R <- response_matrix(adherence_sets(treatments), instrument)
#' KB <- projection_matrices(R, treatments)
#' b <- solve_b_pairs(KB)
#' d <- data.frame(
#'   z = rep(1:2, each = 200),
#'   trt = c(sample(c("a", "b"), 200, TRUE, c(0.7, 0.3)),
#'           sample(c("a", "b"), 200, TRUE, c(0.3, 0.7))),
#'   y = rbinom(400, 1, 0.5)
#' )
#' bootstrap_late(20, d, "z", "trt", "y", KB, b, n_cores = 2)$ci
#' }
#'
#' @export
bootstrap_late <- function(n, data, instrument, treatment, outcome,
                           projections, b, alpha = 0.05,
                           n_cores = default_cores(),
                           cap = FALSE, covariates = NULL, parametric = FALSE,
                           family = NULL, seed = 1234,
                           instrument_levels = NULL, treatments = NULL) {
  datasets <- if (is.data.frame(data)) list(data) else data
  if (!length(datasets) || !all(vapply(datasets, is.data.frame, logical(1)))) {
    stop("`data` must be a data frame or a list of data frames.", call. = FALSE)
  }
  if (is.null(treatments)) treatments <- names(projections)
  if (is.null(instrument_levels)) {
    instrument_levels <- rownames(projections[[1]]$B)
  }

  estimate_once <- function(d) {
    p_z <- estimate_p_z(d, instrument, treatment, covariates, parametric,
                        instrument_levels, treatments)
    q_z <- estimate_q_z(d, instrument, treatment, outcome, covariates,
                        parametric, family, p_z, instrument_levels, treatments)
    unlist(estimate_late(q_z, projections, b, estimate_p_sigma(p_z, projections, b)))
  }

  # Estimate once on the data as given. This surfaces a misspecified call
  # before any workers are started, and fixes the shape and names that every
  # replicate must return, so that a failed replicate can be recorded as NA
  # instead of producing a ragged result.
  template <- estimate_once(datasets[[1]])
  failed <- stats::setNames(rep(NA_real_, length(template)), names(template))

  replicate_once <- function(i, d) {
    resample <- d[sample(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
    tryCatch(estimate_once(resample), error = function(e) failed)
  }

  cluster <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  # Workers start with the default library path, so catIV would be invisible to
  # them whenever it is installed somewhere else, as with renv or a personal
  # library. Hand them the paths this session is using.
  parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
  doParallel::registerDoParallel(cluster)

  # One task per replicate, rather than one per worker. doRNG gives each
  # foreach iteration its own random number stream, so a replicate draws the
  # same resample whatever `n_cores` is. Splitting the work into one chunk per
  # worker instead would tie the streams to the number of workers, and the
  # estimates would depend on the machine they were computed on.
  jobs <- expand.grid(replicate = seq_len(n), dataset = seq_along(datasets))
  job <- NULL  # bound by foreach
  draws <- foreach::foreach(
    job = seq_len(nrow(jobs)), .combine = "cbind", .packages = "catIV",
    .options.RNG = seed
  ) %dorng% {
    replicate_once(jobs$replicate[job], datasets[[jobs$dataset[job]]])
  }
  draws <- matrix(draws, nrow = length(template),
                  dimnames = list(names(template), NULL))

  replicates <- as.data.frame(t(draws))
  effect_range <- range(unlist(lapply(datasets, function(d) d[[outcome]])),
                        na.rm = TRUE)
  cap_value <- diff(effect_range)

  bounds <- lapply(replicates, function(x) {
    keep <- !is.na(x)
    if (cap) keep <- keep & abs(x) <= cap_value
    if (!any(keep)) return(c(NA_real_, NA_real_, 0))
    c(stats::quantile(x[keep], c(alpha / 2, 1 - alpha / 2)), sum(keep))
  })
  ci <- as.data.frame(do.call(rbind, bounds))
  colnames(ci) <- c("lower", "upper", "n_kept")
  rownames(ci) <- colnames(replicates)

  list(alpha = alpha, replicates = replicates, ci = ci)
}

#' A naive instrumental variable estimate ignoring the other treatments
#'
#' Collapses the ordinal instrument to a binary one, indicating whether the
#' first treatment is more encouraged than the second, restricts the data to
#' observations receiving one of the two, and applies the usual Wald ratio.
#'
#' @param contrast A pair of treatments as a string, `"t1_t2"`.
#' @param data A data frame.
#' @param instrument_values The named list of instrument values, as passed to
#'   [response_matrix()].
#' @param instrument,treatment,outcome Column names, as strings.
#'
#' @return A single number: the Wald estimate.
#'
#' @details
#' This is the comparison estimator in paper 1, showing what happens when the
#' remaining treatment alternatives are ignored rather than modelled. It is
#' generally biased under categorical monotonicity, because switching between
#' the two treatments of interest is not the only response to the instrument.
#'
#' Called `NaiveIV()` in the paper code.
#'
#' @export
naive_iv <- function(contrast, data, instrument_values, instrument, treatment,
                     outcome) {
  ts <- split_contrast(contrast)
  encourages_first <- vapply(
    instrument_values,
    function(z) which(z == ts[1]) < which(z == ts[2]),
    logical(1)
  )
  keep <- as.character(data[[treatment]]) %in% ts
  d <- data[keep, , drop = FALSE]
  idx <- match(as.character(d[[instrument]]), names(instrument_values))
  if (anyNA(idx)) {
    stop("Some instrument values in `data` are not present in ",
         "`instrument_values`.", call. = FALSE)
  }
  first <- encourages_first[idx]
  (mean(d[[outcome]][first]) - mean(d[[outcome]][!first])) /
    (mean(as.character(d[[treatment]])[first] == ts[1]) -
       mean(as.character(d[[treatment]])[!first] == ts[1]))
}
