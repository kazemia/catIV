# Report the instrument by name so that the Z column lines up with the row
# names of response_matrix(). Falls back to the index for an unnamed list.
instrument_label <- function(instrument, index) {
  if (is.null(names(instrument))) index else names(instrument)[index]
}

#' Simulate a binary treatment with a binary instrument
#'
#' The two-treatment data-generating process from paper 1: a single binary
#' confounder affects both the treatment and the outcome, and a binary
#' instrument shifts the treatment only.
#'
#' @param n Number of observations.
#' @param z_on_t Effect of the instrument on the treatment.
#' @param v_on_t,v_on_y Effects of the confounder on the treatment and on the
#'   outcome.
#' @param t_on_y Effect of the treatment on the outcome.
#' @param logit If `FALSE` (the default) the effects are probabilities and
#'   combine additively, so the caller must keep every total in `[0, 1]`. If
#'   `TRUE` they are log odds and pass through a logistic link, which is safe
#'   for any values.
#'
#' @return A data frame with columns `Z`, `T` and `Y`, all binary.
#'
#' @details Called `BinarySimulator()` in the paper code, where `logit` was
#'   the argument `OR`.
#'
#' @examples
#' set.seed(1)
#' head(simulate_binary(10, z_on_t = 0.4, v_on_t = 0.3,
#'                      v_on_y = 0.25, t_on_y = 0.2))
#'
#' @export
simulate_binary <- function(n, z_on_t, v_on_t, v_on_y, t_on_y, logit = FALSE) {
  v <- stats::rbinom(n, 1, 0.5)
  z <- stats::rbinom(n, 1, 0.5)
  link <- if (logit) stats::plogis else identity
  t <- stats::rbinom(n, 1, link(v_on_t * v + z_on_t * z))
  y <- stats::rbinom(n, 1, link(v_on_y * v + t_on_y * t))
  data.frame(Z = z, T = t, Y = y)
}

#' Simulate from the adherence set model
#'
#' Generates data the way the identification results assume it arises: a
#' confounder determines the decision team's adherence set, the instrument is
#' drawn independently, and categorical monotonicity then picks the treatment
#' deterministically. Effects estimated from this process are exactly the
#' quantities the package identifies, so it is the right generator for checking
#' that an estimator recovers a known truth.
#'
#' @param n Number of observations.
#' @param instrument A named list of instrument values, as for
#'   [response_matrix()]. Drawn uniformly.
#' @param sets A named list of adherence sets, as for [response_matrix()].
#' @param treatments Character vector of the treatment alternatives.
#' @param v_levels Values the confounder can take.
#' @param v_probs Probability of each confounder level; must sum to one.
#' @param v_to_set A list with one element per confounder level, named by that
#'   level, each a vector of probabilities over `sets` summing to one.
#' @param v_on_y Effect of each confounder level on the outcome probability,
#'   in the order of `v_levels`.
#' @param t_on_y Effect of each treatment on the outcome probability, in the
#'   order of `treatments`.
#' @param choice The choice function, by default [choose_treatment()].
#'
#' @return A data frame with columns `Z` (the name of the instrument value, or
#'   its index if `instrument` is unnamed), `T`, `Y`, and `A` naming the
#'   adherence set each observation was drawn with. `A` is unobservable in
#'   practice and is returned so that simulations can check against it.
#'
#' @details
#' The outcome is Bernoulli with probability `t_on_y[T] + v_on_y[V]`, so the
#' caller must choose effects keeping every total in `[0, 1]`.
#'
#' Called `CatSimulator()` in paper 1's code.
#'
#' @export
simulate_adherence <- function(n, instrument, sets, treatments, v_levels,
                               v_probs, v_to_set, v_on_y, t_on_y,
                               choice = choose_treatment) {
  v <- sample(v_levels, n, replace = TRUE, prob = v_probs)
  set_index <- unlist(lapply(v_to_set[as.character(v)], function(p) {
    sample.int(length(sets), size = 1, prob = p)
  }))
  z_index <- sample.int(length(instrument), size = n, replace = TRUE)
  t <- vapply(seq_len(n),
              function(i) choice(instrument[[z_index[i]]], sets[[set_index[i]]]),
              character(1))
  p_y <- t_on_y[match(t, treatments)] + v_on_y[match(v, v_levels)]
  data.frame(Z = instrument_label(instrument, z_index), T = t,
             Y = stats::rbinom(n, 1, p_y), A = names(sets)[set_index],
             stringsAsFactors = FALSE)
}

#' Simulate a categorical treatment chosen by a multinomial model
#'
#' The generator used for the simulation studies in papers 1 and 2. Unlike
#' [simulate_adherence()], the treatment is drawn from a multinomial model
#' given the instrument and several binary confounders, rather than through a
#' deterministic choice function. Categorical monotonicity therefore holds only
#' approximately, which is the point: it tests the estimators under a process
#' that does not exactly satisfy their assumptions.
#'
#' @param n Number of observations.
#' @param instrument A named list of instrument values, as for
#'   [response_matrix()]. Drawn uniformly.
#' @param treatments Character vector of the treatment alternatives.
#' @param v_on_t A matrix with one row per treatment and one column per
#'   confounder, giving the effect of each confounder on each treatment.
#' @param v_probs Prevalence of each confounder, one value per confounder.
#' @param v_on_y Effect of each confounder on the outcome, one value per
#'   confounder.
#' @param t_on_y Effect of each treatment on the outcome, in the order of
#'   `treatments`.
#' @param intercept Added to every outcome probability.
#' @param observed Logical vector saying which confounders are returned. A
#'   confounder set to `FALSE` acts as unobserved confounding. Recycled to the
#'   number of confounders, so the default returns all of them.
#' @param z_on_t `"exp"` to exponentiate the instrument's encouragement rank
#'   before combining, `"linear"` to use it directly.
#'
#' @return A data frame with columns `Z` (the name of the instrument value, or
#'   its index if `instrument` is unnamed), `T`, `Y`, and one column `V1`,
#'   `V2`, ... per observed confounder.
#'
#' @details
#' Outcome probabilities are clipped to `[0, 1]` after adding the treatment
#' effect, the confounder effects and the intercept.
#'
#' Called `CatSimulator2()` in paper 1's code and `CatSimulator()` in paper
#' 2's; paper 2's version is the same generator with `intercept = 0`.
#'
#' @export
simulate_categorical <- function(n, instrument, treatments, v_on_t, v_probs,
                                 v_on_y, t_on_y, intercept = 0,
                                 observed = TRUE, z_on_t = c("exp", "linear")) {
  z_on_t <- match.arg(z_on_t)
  n_t <- length(treatments)
  n_v <- length(v_probs)
  observed <- rep_len(observed, n_v)

  z_index <- sample.int(length(instrument), size = n, replace = TRUE)
  z_value <- instrument[z_index]

  v <- matrix(NA_real_, nrow = n, ncol = n_v)
  for (j in seq_len(n_v)) v[, j] <- stats::rbinom(n, 1, v_probs[j])

  # Encouragement rank of each treatment under each drawn instrument value.
  rank <- vapply(treatments,
                 function(t) vapply(z_value,
                                    function(z) n_t - which(z == t), numeric(1)),
                 numeric(n))
  rank <- matrix(rank, nrow = n, ncol = n_t,
                 dimnames = list(NULL, treatments))
  if (z_on_t == "exp") rank <- exp(rank)

  p_t <- rank + v %*% t(v_on_t)
  p_t <- p_t / rowSums(p_t)
  t <- apply(p_t, 1, function(p) sample(treatments, size = 1, prob = p))

  p_y <- t_on_y[match(t, treatments)] + as.numeric(v %*% v_on_y) + intercept
  p_y <- pmin(pmax(p_y, 0), 1)

  out <- data.frame(Z = instrument_label(instrument, z_index), T = t,
                    Y = stats::rbinom(n, 1, p_y), stringsAsFactors = FALSE)
  for (j in seq_len(n_v)) {
    if (observed[j]) out[[paste0("V", j)]] <- v[, j]
  }
  out
}
