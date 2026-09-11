#' Heterogeneous instrumental variable estimator of conditional responses
#'
#' The estimator of paper 3. A confounded outcome model is efficient but biased
#' by unobserved confounding; the categorical instrumental variable estimator
#' is unbiased but unstable wherever the target population is small. This fits
#' a weighted model for the difference between the two, and subtracts the
#' fitted bias from the confounded estimate, giving a conditional average
#' treatment response that keeps most of the efficiency of the confounded model
#' without its bias.
#'
#' @param x_values Numeric vector of covariate values to estimate at.
#' @param numerator_model,denominator_model Models of the outcome and of the
#'   treatment, as for [conditional_latr()].
#' @param confounded_model An outcome model adjusting only for observed
#'   covariates, as for [confounded_response()].
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_treatments()].
#' @param density Estimated density of the covariate at `x_values`, as a
#'   numeric vector of the same length. See details.
#' @param instrument_levels,treatments Character vectors of the instrument
#'   values and treatment alternatives.
#' @param instrument,treatment,covariate Column names as the models were
#'   fitted with them.
#' @param density_power Exponent on the density in the regression weights.
#'   Defaults to 1, matching the variance argument in the paper. Use 2 to
#'   reproduce the published analysis code. See the Weights section.
#'
#' @return A data frame with one row per covariate value and treatment:
#'   \describe{
#'     \item{`x`, `treatment`}{where the estimate applies.}
#'     \item{`confounded`}{the covariate-adjusted response.}
#'     \item{`bias`}{the fitted bias.}
#'     \item{`response`}{the debiased conditional average treatment response,
#'       `confounded - bias`.}
#'   }
#'   The fitted bias model is attached as the `"bias_model"` attribute and the
#'   assembled per-solution data as `"bias_data"`.
#'
#' @section Weights:
#' Each observation of the bias enters the weighted regression with weight
#' `density^density_power * p_sigma^2`, where `p_sigma` is the conditional
#' probability of belonging to the target population, that is the denominator
#' of the instrumental variable estimator.
#'
#' The paper specifies the weight as the square of that denominator times the
#' density of the covariate, consistent with its variance expression, in which
#' the variance is proportional to `1 / (p_sigma^2 * density)`. That gives
#' `density_power = 1`, which is the default here.
#'
#' The analysis code used for the published clinical results computed the
#' weight with the density **squared**. Set `density_power = 2` to reproduce
#' those numbers exactly. The default does not reproduce them, because the
#' squared density was not what the paper specifies.
#'
#' @section Contrasts:
#' Responses are estimated one treatment at a time. Form a conditional average
#' treatment effect by differencing two of them, which [hiv_contrast()] does.
#'
#' @details
#' This estimator is assembled from the analysis script for paper 3 rather than
#' from its function file. Its components, [conditional_p_sigma()] and
#' [conditional_latr()], reproduce the originals exactly.
#'
#' Inference is by bootstrapping the whole procedure; there is no closed form.
#'
#' @export
hiv_cate <- function(x_values, numerator_model, denominator_model,
                     confounded_model, projections, b, density,
                     instrument_levels, treatments, instrument = "Z",
                     treatment = "T", covariate = "X", density_power = 1) {
  if (length(density) != length(x_values)) {
    stop("`density` must have one value per element of `x_values`.",
         call. = FALSE)
  }

  p_sigma <- lapply(x_values, conditional_p_sigma, denominator_model,
                    projections, b, instrument_levels, treatments, instrument,
                    covariate)
  latr <- lapply(x_values, conditional_latr, numerator_model,
                 denominator_model, projections, b, instrument_levels,
                 treatments, instrument, treatment, covariate)
  confounded <- confounded_response(x_values, confounded_model, treatments,
                                    treatment, covariate)
  conf_at <- function(t) {
    confounded$response[match(paste(x_values, t),
                              paste(confounded$x, confounded$treatment))]
  }

  blocks <- list()
  for (t in treatments) {
    for (s in rownames(b[[t]])) {
      pi_hat <- vapply(p_sigma, function(v) v[[t]][[s]], numeric(1))
      e_hat <- vapply(latr, function(v) v[[t]][[s]], numeric(1))
      c_hat <- conf_at(t)
      blocks[[length(blocks) + 1L]] <- data.frame(
        x = x_values, treatment = t, solution = s,
        p_sigma = pi_hat, latr = e_hat, confounded = c_hat,
        density = density, bias = c_hat - e_hat,
        w = density^density_power * pi_hat^2,
        stringsAsFactors = FALSE
      )
    }
  }
  bias_data <- do.call(rbind, blocks)

  bias_model <- stats::lm(bias ~ treatment * x + x * solution,
                          weights = bias_data$w, data = bias_data)

  solutions <- unique(bias_data$solution)
  out <- lapply(treatments, function(t) {
    weighted <- lapply(solutions, function(s) {
      fitted_bias <- stats::predict(
        bias_model,
        newdata = data.frame(x = x_values, treatment = t, solution = s,
                             stringsAsFactors = FALSE)
      )
      sub <- bias_data[bias_data$solution == s, , drop = FALSE]
      # A solution can be identified for more than one treatment, in which case
      # its weights are averaged across them at each covariate value.
      wts <- if (nrow(sub) > length(x_values)) {
        as.numeric(tapply(sqrt(sub$w), factor(sub$x, levels = x_values), mean))
      } else {
        sqrt(sub$w)
      }
      list(numerator = wts * fitted_bias, denominator = wts)
    })
    bias_hat <- Reduce(`+`, lapply(weighted, `[[`, "numerator")) /
      Reduce(`+`, lapply(weighted, `[[`, "denominator"))
    c_hat <- conf_at(t)
    data.frame(x = x_values, treatment = t, confounded = c_hat,
               bias = bias_hat, response = c_hat - bias_hat,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  attr(out, "bias_model") <- bias_model
  attr(out, "bias_data") <- bias_data
  out
}

#' Contrast two treatments from a heterogeneous IV fit
#'
#' @param fit Output of [hiv_cate()].
#' @param t1,t2 The two treatments to contrast, as strings. The result is the
#'   response under `t1` minus the response under `t2`.
#'
#' @return A data frame with columns `x`, `contrast` and `estimate`.
#'
#' @export
hiv_contrast <- function(fit, t1, t2) {
  a <- fit[fit$treatment == t1, , drop = FALSE]
  b <- fit[fit$treatment == t2, , drop = FALSE]
  if (!nrow(a) || !nrow(b)) {
    stop("Both treatments must be present in `fit`.", call. = FALSE)
  }
  b <- b[match(a$x, b$x), , drop = FALSE]
  data.frame(x = a$x, contrast = paste(t1, "-", t2),
             estimate = a$response - b$response,
             stringsAsFactors = FALSE)
}

#' Average a conditional effect over the covariate distribution
#'
#' Collapses a conditional average treatment effect to a single number by
#' averaging over the covariate, weighted by its estimated density.
#'
#' @param contrast Output of [hiv_contrast()].
#' @param density Estimated density of the covariate at `contrast$x`.
#'
#' @return A single number.
#'
#' @export
average_contrast <- function(contrast, density) {
  if (length(density) != nrow(contrast)) {
    stop("`density` must have one value per row of `contrast`.", call. = FALSE)
  }
  sum(contrast$estimate * density) / sum(density)
}
