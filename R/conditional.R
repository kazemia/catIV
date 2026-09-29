# Build the newdata grid used to read a fitted model at one covariate value.
conditional_grid <- function(x, instrument_levels, instrument, covariate) {
  grid <- data.frame(
    factor(instrument_levels, levels = instrument_levels),
    rep(x, length(instrument_levels))
  )
  names(grid) <- c(instrument, covariate)
  grid
}

#' Conditional treatment probabilities from a fitted model
#'
#' Reads `P(T = t | Z = z, X = x)` off a model of treatment given the
#' instrument and an effect-modifying covariate. This is the conditional
#' counterpart of [estimate_p_z()], and the denominator of the categorical
#' instrumental variable estimator.
#'
#' @param x A single value of the effect-modifying covariate.
#' @param model A fitted model of treatment on the instrument and the
#'   covariate, such as [nnet::multinom()] or a random forest. Anything whose
#'   `predict()` method returns class probabilities will work.
#' @param instrument_levels Character vector of the instrument values, in the
#'   order used to build the projection matrices.
#' @param treatments Character vector of the treatment alternatives.
#' @param instrument,covariate Names of the instrument and covariate columns as
#'   the model was fitted with them.
#'
#' @return A numeric matrix with one row per instrument value and one column
#'   per treatment, the same shape as [estimate_p_z()] returns.
#'
#' @details
#' A two-class model whose `predict()` returns a single vector of
#' probabilities is expanded to two columns, taking the vector to be the
#' probability of the second treatment.
#'
#' @export
conditional_p_z <- function(x, model, instrument_levels, treatments,
                            instrument = "Z", covariate = "X") {
  grid <- conditional_grid(x, instrument_levels, instrument, covariate)
  type <- if (inherits(model, "randomForest")) "prob" else "probs"
  probs <- stats::predict(model, newdata = grid, type = type)

  if (is.null(dim(probs)) || ncol(as.matrix(probs)) == 1L) {
    # Two-class model: the vector holds the probability of the second level.
    second <- as.numeric(probs)
    probs <- cbind(1 - second, second)
    colnames(probs) <- model_levels(model, treatments)
  }
  probs <- as.matrix(probs)

  out <- matrix(0, nrow = length(instrument_levels), ncol = length(treatments),
                dimnames = list(instrument_levels, treatments))
  shared <- intersect(colnames(probs), treatments)
  out[, shared] <- probs[, shared, drop = FALSE]
  out
}

model_levels <- function(model, treatments) {
  lev <- tryCatch(model$lev, error = function(e) NULL)
  if (is.null(lev) || length(lev) != 2L) treatments[1:2] else lev
}

#' Conditional probability of belonging to each target population
#'
#' The conditional counterpart of [estimate_p_sigma_treatments()]: for each
#' treatment and each solution `b`, the probability
#' `P(b[A] = 1 | X = x)` of a decision team being in the identified target
#' population at covariate value `x`.
#'
#' @inheritParams conditional_p_z
#' @param projections Output of [projection_matrices()].
#' @param b Output of [solve_b_treatments()].
#'
#' @return A named list over treatments of named numeric vectors over the
#'   solutions in `b`.
#'
#' @details
#' Positivity requires this to be strictly positive for the covariate values of
#' interest. Where it approaches zero the instrumental variable estimator has a
#' vanishing denominator and becomes unstable, which is the problem the
#' estimator in [hiv_cate()] is designed to work around.
#'
#' Called `Pi_prob_estimator()` in the paper 3 code.
#'
#' @export
conditional_p_sigma <- function(x, model, projections, b, instrument_levels,
                                treatments, instrument = "Z",
                                covariate = "X") {
  p_z <- conditional_p_z(x, model, instrument_levels, treatments, instrument,
                         covariate)
  estimate_p_sigma_treatments(p_z, projections, b)
}

#' Conditional local average treatment responses
#'
#' The categorical instrumental variable estimator of paper 3: for each
#' treatment and each solution `b`, estimates
#' `E(Y^{T=t} | b[A] = 1, X = x)`.
#'
#' @inheritParams conditional_p_sigma
#' @param numerator_model A fitted model of the outcome on the instrument, the
#'   treatment and the covariate, for example
#'   `lm(y ~ Z * trt * X, data = d)`.
#' @param denominator_model A fitted model of the treatment on the instrument
#'   and the covariate, as for [conditional_p_z()].
#' @param treatment Name of the treatment column as the numerator model was
#'   fitted with it.
#'
#' @return A named list over treatments of named numeric vectors over the
#'   solutions in `b`.
#'
#' @details
#' The numerator `Q_Z(t, X)` is formed by multiplying the modelled outcome
#' means `E(Y | Z, T, X)` by the modelled treatment probabilities
#' `P(T | Z, X)`, then projecting through `b B_t^+`, exactly as in the
#' unconditional case. Neither model is a local estimator, so the variance
#' expression used to derive the weights in [hiv_cate()] is an approximation
#' in this setting.
#'
#' Called `CIV_estimator()` in the paper 3 code.
#'
#' @export
conditional_latr <- function(x, numerator_model, denominator_model,
                             projections, b, instrument_levels, treatments,
                             instrument = "Z", treatment = "T",
                             covariate = "X") {
  p_z <- conditional_p_z(x, denominator_model, instrument_levels, treatments,
                         instrument, covariate)

  grid_args <- list(model = numerator_model, x)
  names(grid_args)[2] <- covariate
  grid_args[[instrument]] <- factor(instrument_levels,
                                    levels = instrument_levels)
  grid_args[[treatment]] <- treatments
  grid <- do.call(marginaleffects::datagrid, grid_args)
  long <- marginaleffects::predictions(numerator_model, newdata = grid)

  means <- widen_effects(long, instrument, treatment, "estimate",
                         instrument_levels, treatments)
  means[is.na(means)] <- 0
  q_z <- means * p_z

  p_sigma <- estimate_p_sigma_treatments(p_z, projections, b)
  estimate_latr(q_z, projections, b, p_sigma)
}

#' Conditional treatment responses from a confounded outcome model
#'
#' Reads `E(Y | T = t, X = x)` off an outcome model that adjusts only for
#' observed covariates. These estimates are efficient but biased by any
#' unobserved confounding; [hiv_cate()] combines them with the instrumental
#' variable estimates to remove that bias.
#'
#' @param x Values of the effect-modifying covariate.
#' @param model A fitted outcome model, for example a `glm` or a `mgcv::gam`.
#' @param treatments Character vector of the treatment alternatives.
#' @param treatment,covariate Names of the treatment and covariate columns as
#'   the model was fitted with them.
#'
#' @return A data frame with columns `x`, `treatment` and `response`.
#'
#' @details
#' Covariates other than `x` are held at typical values, numeric ones at their
#' mean and categorical ones at their mode, as `marginaleffects::datagrid()`
#' does. They are not averaged over their distribution.
#'
#' Called `Confounded_estimator()` and `CatConfounded_estimator()` in the paper
#' code.
#'
#' @export
confounded_response <- function(x, model, treatments, treatment = "T",
                                covariate = "X") {
  grid_args <- list(model = model, x)
  names(grid_args)[2] <- covariate
  grid_args[[treatment]] <- factor(treatments, levels = treatments)
  grid <- do.call(marginaleffects::datagrid, grid_args)
  long <- as.data.frame(
    marginaleffects::avg_predictions(model, newdata = grid,
                                     by = c(covariate, treatment))
  )
  data.frame(
    x = as.numeric(as.character(long[[covariate]])),
    treatment = as.character(long[[treatment]]),
    response = long$estimate,
    stringsAsFactors = FALSE
  )
}
