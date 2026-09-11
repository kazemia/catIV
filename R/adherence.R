#' Enumerate all possible adherence sets
#'
#' An adherence set is the set of treatment alternatives that a decision team
#' would consider for a given patient, that is the alternatives with non-zero
#' probability of being chosen once the effect of the instrument is
#' marginalised out. With `NT` treatment alternatives there are `2^NT - 1`
#' non-empty adherence sets.
#'
#' @param treatments Character vector of the treatment alternatives,
#'   `Supp(T)`.
#'
#' @return A named list of character vectors, one per adherence set, named
#'   `A1`, `A2`, ... in order of increasing set size.
#'
#' @details
#' In applications the full enumeration is often implausible, and a smaller
#' list of clinically realistic adherence sets is supplied to
#' [response_matrix()] instead. Restricting the list makes more effects
#' identifiable, at the cost of an additional substantive assumption.
#'
#' @examples
#' adherence_sets(c("a", "b", "c"))
#'
#' @export
adherence_sets <- function(treatments) {
  treatments <- as.character(treatments)
  if (anyDuplicated(treatments)) {
    stop("`treatments` must not contain duplicates.", call. = FALSE)
  }
  if (length(treatments) < 2L) {
    stop("`treatments` must contain at least two alternatives.", call. = FALSE)
  }
  sets <- unlist(
    lapply(seq_along(treatments), utils::combn, x = treatments,
           simplify = FALSE),
    recursive = FALSE
  )
  names(sets) <- paste0("A", seq_along(sets))
  sets
}

#' The choice function implied by categorical monotonicity
#'
#' Categorical monotonicity states that the treatment received is a
#' deterministic function of the instrument and the adherence set: the decision
#' team picks the alternative in the adherence set that the instrument
#' encourages most.
#'
#' @param z Character vector of all treatment alternatives, sorted in
#'   *decreasing* order of encouragement by the instrument. In the pricing
#'   application this is the alternatives ordered from cheapest to most
#'   expensive.
#' @param a Character vector of the alternatives the decision team adheres to,
#'   i.e. one adherence set.
#'
#' @return A length-one character vector: the element of `a` appearing earliest
#'   in `z`.
#'
#' @examples
#' choose_treatment(c("a", "b", "c"), c("b", "c"))
#'
#' @export
choose_treatment <- function(z, a) {
  position <- match(a, z)
  if (all(is.na(position))) {
    stop("No element of the adherence set appears in the instrument value.",
         call. = FALSE)
  }
  a[which.min(position)]
}

#' Build the response matrix
#'
#' The response matrix records which treatment each adherence set receives
#' under each value of the instrument. It is the design matrix from which the
#' indicator matrices `B_t` are formed.
#'
#' @param sets A named list of adherence sets, typically from
#'   [adherence_sets()] or a hand-specified list of clinically realistic sets.
#' @param instrument A named list of instrument values. Each element is a
#'   character vector of all treatment alternatives sorted in decreasing order
#'   of encouragement.
#' @param choice The choice function. Defaults to [choose_treatment()], the
#'   function implied by categorical monotonicity.
#'
#' @return A data frame with one row per instrument value and one column per
#'   adherence set, holding the treatment received.
#'
#' @examples
#' instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"))
#' response_matrix(adherence_sets(c("a", "b", "c")), instrument)
#'
#' @export
response_matrix <- function(sets, instrument, choice = choose_treatment) {
  if (!length(sets) || !length(instrument)) {
    stop("`sets` and `instrument` must both be non-empty.", call. = FALSE)
  }
  out <- vapply(
    sets,
    function(a) vapply(instrument, function(z) choice(z, a), character(1L)),
    character(length(instrument))
  )
  out <- matrix(out, nrow = length(instrument), ncol = length(sets),
                dimnames = list(names(instrument), names(sets)))
  as.data.frame(out)
}
