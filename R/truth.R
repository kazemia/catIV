#' True principal stratum probabilities for a binary instrument
#'
#' Computes the exact conditional probabilities of each principal stratum, and
#' of receiving treatment, for the binary data-generating process used in the
#' simulations of paper 3. Because these are the values the estimators are
#' trying to recover, they give a simulation something to be checked against.
#'
#' The process is
#' `P(T = 1 | Z, U, V, X) = plogis(z_on_t * Z + u_on_t * U + v_on_t * V + x_on_t * X)`
#' with `U`, `V` and `Z` independent Bernoulli variables. `U` is unobserved
#' confounding, `V` is observed confounding.
#'
#' @param x Values of the effect-modifying covariate.
#' @param z_on_t,x_on_t,u_on_t,v_on_t Log-odds effects of the instrument, the
#'   covariate, the unobserved confounder and the observed confounder on the
#'   treatment.
#' @param u_prob,v_prob,z_prob Prevalences of the unobserved confounder, the
#'   observed confounder and the instrument.
#'
#' @return A data frame with one row per value of `x`:
#'   \describe{
#'     \item{`complier`}{`P(T = 1 | Z = 1, X) - P(T = 1 | Z = 0, X)`.}
#'     \item{`always_taker`}{`P(T = 1 | Z = 0, X)`.}
#'     \item{`never_taker`}{`P(T = 0 | Z = 1, X)`.}
#'     \item{`treated`}{`P(T = 1 | X)`, marginal over the instrument.}
#'   }
#'   The three stratum probabilities sum to one.
#'
#' @details
#' This replaces `C_probability()`, `AT_probability()`, `NT_probability()` and
#' `T_probability()` from the paper 3 code. Those were written as four separate
#' expressions, and `T_probability()` took no `VT` or `VP` arguments, so it
#' marginalised over the unobserved confounder and the instrument but not over
#' the observed confounder. Setting `v_on_t = 0` here reproduces that
#' behaviour; any other value corrects it. Everything is marginalised over the
#' same 2x2 grid of `U` and `V`, so the four quantities cannot drift apart.
#'
#' @examples
#' binary_iv_truth(c(0, 1, 2), z_on_t = 1.5, x_on_t = 0.4, u_on_t = 1,
#'                 v_on_t = 0.8, u_prob = 0.5, v_prob = 0.5, z_prob = 0.5)
#'
#' @export
binary_iv_truth <- function(x, z_on_t, x_on_t, u_on_t, v_on_t, u_prob, v_prob,
                            z_prob = 0.5) {
  # P(T = 1 | Z = z, X = x), marginal over U and V.
  treated_given_z <- function(z) {
    total <- 0
    for (u in c(0, 1)) {
      for (v in c(0, 1)) {
        share <- (if (u == 1) u_prob else 1 - u_prob) *
          (if (v == 1) v_prob else 1 - v_prob)
        total <- total +
          share * stats::plogis(z_on_t * z + u_on_t * u + v_on_t * v +
                                  x_on_t * x)
      }
    }
    total
  }
  encouraged <- treated_given_z(1)
  not_encouraged <- treated_given_z(0)
  data.frame(
    x = x,
    complier = encouraged - not_encouraged,
    always_taker = not_encouraged,
    never_taker = 1 - encouraged,
    treated = z_prob * encouraged + (1 - z_prob) * not_encouraged
  )
}
