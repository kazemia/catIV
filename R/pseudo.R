#' Build the pseudo-population behind an identified effect
#'
#' Attaches a weight to each observation receiving one of the two treatments
#' in a contrast, such that a weighted comparison of outcomes reproduces the
#' estimate from [estimate_late()]. The weighted sample is the target trial
#' the instrumental variable analysis emulates, so its weighted baseline
#' characteristics describe the population the effect applies to.
#'
#' @param contrast A pair of treatments as a string, `"t1_t2"`.
#' @param data A data frame.
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_pairs()].
#' @param p_sigma Output of [estimate_p_sigma()].
#' @param instrument,treatment Column names, as strings.
#' @param solution Which identified target population to build the
#'   pseudo-population for, given as a solution name such as `"A4+A7"`, or as
#'   a positive integer indexing the rows of `b[[contrast]]`. Defaults to the
#'   first.
#'
#' @return The subset of `data` receiving one of the two treatments, with an
#'   added numeric column `w` holding the weights.
#'
#' @details
#' `p_sigma` must be the unadjusted estimate, computed from the same
#' `estimate_p_z()` output as the data being weighted. The weights are built
#' from raw cell proportions, so passing a covariate-standardised `P_Sigma`
#' breaks the correspondence: the weighted comparison no longer reproduces the
#' estimate exactly and the weights no longer average to one.
#'
#' Given the unadjusted `p_sigma`, the weights average to exactly one within
#' each treatment arm. Both `weighted.mean(x, w)` and `mean(x * w)` therefore
#' give the same weighted average of a covariate within an arm, which is what
#' makes the weighted baseline table meaningful.
#'
#' Weights can be negative. `b B_t^+` is a contrast across instrument values
#' and need not be non-negative, so some observations enter with a negative
#' sign. This is expected, not a numerical failure.
#'
#' Called `PseudoPopulator()` in the paper code, where the target population
#' was selected by the integer `Pi_index`.
#'
#' @export
pseudo_population <- function(contrast, data, projections, b, p_sigma,
                              instrument, treatment, solution = 1L) {
  ts <- split_contrast(contrast)
  b_contrast <- b[[contrast]]
  if (is.null(b_contrast) || nrow(b_contrast) == 0L) {
    stop("No effect is identified for contrast ", sQuote(contrast), ".",
         call. = FALSE)
  }
  row <- resolve_solution(b_contrast, solution, contrast)

  z_chr <- as.character(data[[instrument]])
  z_levels <- rownames(projections[[1]]$B)
  idx <- match(z_chr, z_levels)
  if (anyNA(idx)) {
    stop("Some instrument values in `data` are not among the instrument ",
         "values used to build the projection matrices.", call. = FALSE)
  }

  # One over the share of the sample at each instrument value.
  share <- as.numeric(table(factor(z_chr, levels = z_levels)))
  w_z <- nrow(data) / share

  keep <- as.character(data[[treatment]]) %in% ts
  pseudo <- data[keep, , drop = FALSE]
  pseudo_idx <- idx[keep]
  pseudo_trt <- as.character(pseudo[[treatment]])

  w <- w_z[pseudo_idx]
  for (k in seq_along(ts)) {
    t_k <- ts[k]
    arm <- pseudo_trt == t_k
    if (!any(arm)) next
    beta <- as.numeric(b_contrast[row, , drop = FALSE] %*%
                         projections[[t_k]]$B_plus)
    w_t <- sum(as.character(data[[treatment]]) == t_k) / nrow(data)
    w_sigma <- 1 / p_sigma[[if (k == 1L) "arm1" else "arm2"]][[contrast]][row]
    w[arm] <- w[arm] * w_t * beta[pseudo_idx[arm]] * w_sigma
  }
  pseudo$w <- w
  pseudo
}

# Accept a solution either by name or by row index.
resolve_solution <- function(b_contrast, solution, contrast) {
  if (is.character(solution)) {
    row <- match(solution, rownames(b_contrast))
    if (is.na(row)) {
      stop("Solution ", sQuote(solution), " is not identified for contrast ",
           sQuote(contrast), ". Available: ",
           paste(sQuote(rownames(b_contrast)), collapse = ", "), ".",
           call. = FALSE)
    }
    return(row)
  }
  row <- as.integer(solution)
  if (is.na(row) || row < 1L || row > nrow(b_contrast)) {
    stop("`solution` must be between 1 and ", nrow(b_contrast),
         " for contrast ", sQuote(contrast), ".", call. = FALSE)
  }
  row
}
