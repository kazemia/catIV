# Split a contrast name such as "Inf_Ada" into its two treatments.
split_contrast <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  if (length(parts) != 2L) {
    stop("Contrast name ", sQuote(x), " does not have the form ",
         sQuote("treatment1_treatment2"), ".", call. = FALSE)
  }
  parts
}

# b %*% B_t^+ %*% m, returned as a plain named vector over the solutions in b.
project <- function(b, projections, treatment, m) {
  out <- as.numeric(b %*% (projections[[treatment]]$B_plus %*% m))
  stats::setNames(out, rownames(b))
}

#' Probability of belonging to each identified target population
#'
#' For every solution `b`, `b B_t^+ P_Z(t)` is the probability that a decision
#' team's adherence set satisfies `b[A] = 1`. Because a contrast involves two
#' treatments, the probability can be computed from either arm; the two agree
#' in the population but not in a finite sample.
#'
#' @param p_z A matrix from [estimate_p_z()].
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_pairs()].
#'
#' @return A list of three elements, each a named list over contrasts of named
#'   numeric vectors over the solutions in `b`:
#'   \describe{
#'     \item{`arm1`}{computed from the first treatment of the pair.}
#'     \item{`arm2`}{computed from the second treatment of the pair.}
#'     \item{`average`}{the weighted average of the two, weighting each arm by
#'       its marginal share `colSums(p_z)`.}
#'   }
#'
#' @details
#' `average` is the quantity used to screen out contrasts whose target
#' population is too small to estimate reliably. Called `P_SigmaIdentifier()`
#' in the paper code, where the three elements were named `E1`, `E2` and `WA`.
#'
#' @export
estimate_p_sigma <- function(p_z, projections, b) {
  weights <- colSums(p_z)
  arm1 <- arm2 <- average <- vector("list", length(b))
  names(arm1) <- names(arm2) <- names(average) <- names(b)
  for (nm in names(b)) {
    ts <- split_contrast(nm)
    a1 <- project(b[[nm]], projections, ts[1], p_z[, ts[1]])
    a2 <- project(b[[nm]], projections, ts[2], p_z[, ts[2]])
    w1 <- weights[[ts[1]]]
    w2 <- weights[[ts[2]]]
    arm1[[nm]] <- a1
    arm2[[nm]] <- a2
    average[[nm]] <- (a1 * w1 + a2 * w2) / (w1 + w2)
  }
  list(arm1 = arm1, arm2 = arm2, average = average)
}

#' Identify local average treatment effects
#'
#' For each contrast and each solution `b`, estimates
#' `E(Y^{T=t} - Y^{T=t'} | b[A] = 1)`, the effect among decision teams whose
#' adherence set satisfies `b[A] = 1`.
#'
#' @param q_z A matrix from [estimate_q_z()].
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_pairs()].
#' @param p_sigma Output of [estimate_p_sigma()].
#' @param scale `"difference"` for a difference in means, `"ratio"` for a risk
#'   ratio.
#' @param denominator Which target-population probability to divide the two
#'   arms by: `"arm"` uses each arm's own estimate, `"average"` uses the
#'   weighted average for both.
#' @param arms If `TRUE`, also return the two local average treatment
#'   *responses* that the contrast is built from.
#'
#' @return With `arms = FALSE` (the default), a named list over contrasts of
#'   named numeric vectors over the solutions in `b`. With `arms = TRUE`, a
#'   list of three such objects: `arm1`, `arm2` and `effect`.
#'
#' @details
#' Note that with `scale = "ratio"` and `denominator = "average"` the shared
#' denominator cancels, so the result is the ratio of the raw projected
#' quantities.
#'
#' Called `LATEIdentifier()` in the paper code, where `scale` and
#' `denominator` were the logical arguments `RR` and `AverageProb`.
#'
#' @export
estimate_late <- function(q_z, projections, b, p_sigma,
                          scale = c("difference", "ratio"),
                          denominator = c("arm", "average"),
                          arms = FALSE) {
  scale <- match.arg(scale)
  denominator <- match.arg(denominator)

  arm1 <- arm2 <- effect <- vector("list", length(b))
  names(arm1) <- names(arm2) <- names(effect) <- names(b)
  for (nm in names(b)) {
    ts <- split_contrast(nm)
    n1 <- project(b[[nm]], projections, ts[1], q_z[, ts[1]])
    n2 <- project(b[[nm]], projections, ts[2], q_z[, ts[2]])
    if (denominator == "arm") {
      d1 <- p_sigma$arm1[[nm]]
      d2 <- p_sigma$arm2[[nm]]
    } else {
      d1 <- d2 <- p_sigma$average[[nm]]
    }
    arm1[[nm]] <- n1 / d1
    arm2[[nm]] <- n2 / d2
    effect[[nm]] <- if (scale == "difference") {
      arm1[[nm]] - arm2[[nm]]
    } else {
      arm1[[nm]] / arm2[[nm]]
    }
  }
  if (!arms) return(effect)
  list(arm1 = arm1, arm2 = arm2, effect = effect)
}

#' Identify local average treatment responses
#'
#' The single-arm counterpart of [estimate_late()]: for each treatment and each
#' solution `b`, estimates `E(Y^{T=t} | b[A] = 1)`.
#'
#' @param q_z A matrix from [estimate_q_z()].
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_treatments()], indexed by single treatments
#'   rather than by contrasts.
#' @param p_sigma A named list over treatments of named numeric vectors, as
#'   returned by [estimate_p_sigma_treatments()].
#'
#' @return A named list over treatments of named numeric vectors over the
#'   solutions in `b`.
#'
#' @details Called `LATOIdentifier()` in the paper code.
#'
#' @export
estimate_latr <- function(q_z, projections, b, p_sigma) {
  out <- vector("list", length(b))
  names(out) <- names(b)
  for (nm in names(b)) {
    out[[nm]] <- project(b[[nm]], projections, nm, q_z[, nm]) / p_sigma[[nm]]
  }
  out
}

#' Probability of belonging to each single-treatment target population
#'
#' The single-arm counterpart of [estimate_p_sigma()], for use with
#' [solve_b_treatments()] and [estimate_latr()].
#'
#' @param p_z A matrix from [estimate_p_z()].
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_treatments()].
#'
#' @return A named list over treatments of named numeric vectors over the
#'   solutions in `b`.
#'
#' @details This is the version of `P_SigmaIdentifier()` defined in
#'   `HTEfunctions.R` for paper 3, which indexes `b` by single treatments.
#'
#' @export
estimate_p_sigma_treatments <- function(p_z, projections, b) {
  out <- vector("list", length(b))
  names(out) <- names(b)
  for (nm in names(b)) {
    out[[nm]] <- project(b[[nm]], projections, nm, p_z[, nm])
  }
  out
}
