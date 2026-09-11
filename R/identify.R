#' Indicator matrices and null-space projectors
#'
#' For each treatment alternative `t` this builds the binary indicator matrix
#' `B_t`, where `B_t[i, n]` is one when adherence set `n` receives treatment
#' `t` under instrument value `i`; its Moore-Penrose pseudo-inverse `B_t^+`;
#' and the projector `K_t = I - B_t^+ B_t` onto the null space of `B_t`.
#'
#' The vectors `b` that identify a local average treatment response are the
#' binary solutions of `b = b B_t^+ B_t`, equivalently `b K_t = 0`. See
#' [solve_b_treatments()] and [solve_b_pairs()].
#'
#' @param response A response matrix from [response_matrix()].
#' @param treatments Character vector of the treatment alternatives.
#' @param tolerance Number of decimal places to which `K_t` is rounded. This
#'   removes floating-point noise from the pseudo-inverse. It must be *larger*
#'   than the `tolerance` later passed to the solvers, so that entries rounded
#'   to zero here are also treated as zero there.
#'
#' @return A named list with one element per treatment, each a list of
#'   `K`, `B` and `B_plus`.
#'
#' @examples
#' instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"),
#'                    z3 = c("b", "c", "a"), z4 = c("c", "b", "a"))
#' treatments <- c("a", "b", "c")
#' R <- response_matrix(adherence_sets(treatments), instrument)
#' KB <- projection_matrices(R, treatments)
#' KB[["a"]]$B
#'
#' @export
projection_matrices <- function(response, treatments, tolerance = 4) {
  response <- as.matrix(response)
  n_sets <- ncol(response)
  out <- lapply(treatments, function(t) {
    B <- 1 * (response == t)
    B_plus <- MASS::ginv(B)
    K <- round(diag(n_sets) - B_plus %*% B, tolerance)
    dimnames(K) <- list(colnames(response), colnames(response))
    dimnames(B_plus) <- list(colnames(response), rownames(response))
    list(K = K, B = B, B_plus = B_plus)
  })
  names(out) <- treatments
  out
}

# Enumerate every binary b satisfying `b K = 0` for all K in `Ks`.
#
# Variables whose column is entirely zero across the stacked constraints are
# unconstrained and are enumerated freely. Constraint rows with a single
# non-zero entry force that variable to zero. Both reductions are exact; the
# remaining variables are enumerated exhaustively and every candidate is
# re-checked against the full stack, so the returned set is complete.
solve_b_stack <- function(Ks, tolerance) {
  K_stack <- do.call(rbind, Ks)
  n_sets <- ncol(K_stack)
  set_names <- colnames(Ks[[1]])
  eps <- 10^(-tolerance)

  nonzero <- K_stack != 0
  free <- which(colSums(nonzero) == 0L)
  single <- which(rowSums(nonzero) == 1L)
  forced_zero <- unique(max.col(nonzero[single, , drop = FALSE], ties.method = "first"))

  enumerated <- setdiff(seq_len(n_sets), union(free, forced_zero))
  n_unknown <- length(enumerated) + length(free)
  if (n_unknown > 24L) {
    stop("The reduced problem has ", n_unknown, " unconstrained binary ",
         "variables, which is too many to enumerate exhaustively. Consider ",
         "restricting the list of adherence sets.", call. = FALSE)
  }

  grid <- as.matrix(expand.grid(rep(list(c(0, 1)), length(enumerated))))
  b <- matrix(0, nrow = max(nrow(grid), 1L), ncol = n_sets)
  if (length(enumerated)) b[, enumerated] <- grid

  # Full cross-product over the unconstrained variables.
  for (fv in free) {
    with_one <- b
    with_one[, fv] <- 1
    b <- rbind(b, with_one)
  }

  # `K_stack` is the vertical stack of the constraints, so post-multiplying by
  # its transpose evaluates b %*% K for every constraint block at once.
  satisfied <- rowSums(abs(tcrossprod(b, K_stack)) > eps) == 0L
  b <- b[satisfied & rowSums(b) > 0, , drop = FALSE]
  b <- unique(b)
  colnames(b) <- set_names
  b <- b[order(rowSums(b), apply(b, 1, paste, collapse = "")), , drop = FALSE]
  # Name each solution by the adherence sets it selects, so that downstream
  # results are addressed by meaning rather than by position.
  rownames(b) <- apply(b, 1L, function(row) {
    paste(set_names[row == 1], collapse = "+")
  })
  b
}

#' Identify the target populations of local average treatment responses
#'
#' Finds every binary vector `b` solving `b = b B_t^+ B_t` for each treatment
#' alternative separately. Each solution defines a target population, namely
#' the decision teams whose adherence set `a` satisfies `b[a] = 1`, for which
#' the local average treatment response `E(Y^{T=t} | b[A] = 1)` is identified.
#'
#' This is the single-arm quantity used by the conditional estimators. For
#' contrasts between two treatments use [solve_b_pairs()].
#'
#' @param projections Output of [projection_matrices()].
#' @param tolerance Number of decimal places of accuracy. Must be *smaller*
#'   than the `tolerance` used in [projection_matrices()].
#'
#' @return A named list with one matrix per treatment. Each row is a solution
#'   `b`, each column an adherence set.
#'
#' @examples
#' instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"),
#'                    z3 = c("b", "c", "a"), z4 = c("c", "b", "a"))
#' treatments <- c("a", "b", "c")
#' R <- response_matrix(adherence_sets(treatments), instrument)
#' solve_b_treatments(projection_matrices(R, treatments))
#'
#' @export
solve_b_treatments <- function(projections, tolerance = 3) {
  out <- lapply(projections, function(x) solve_b_stack(list(x$K), tolerance))
  names(out) <- names(projections)
  out
}

#' Identify the target populations of local average treatment effects
#'
#' Finds every binary vector `b` that solves `b = b B_t^+ B_t = b B_t'^+ B_t'`
#' simultaneously, for each pair of treatment alternatives. Each solution
#' defines a target population for which the contrast
#' `E(Y^{T=t} - Y^{T=t'} | b[A] = 1)` is identified.
#'
#' @param projections Output of [projection_matrices()].
#' @param tolerance Number of decimal places of accuracy. Must be *smaller*
#'   than the `tolerance` used in [projection_matrices()].
#'
#' @return A named list with one matrix per pair of treatments, named
#'   `"t1_t2"`. Each row is a solution `b`, each column an adherence set. A
#'   pair with no solution yields a zero-row matrix, meaning no effect is
#'   identified for that contrast.
#'
#' @examples
#' instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"),
#'                    z3 = c("b", "c", "a"), z4 = c("c", "b", "a"))
#' treatments <- c("a", "b", "c")
#' R <- response_matrix(adherence_sets(treatments), instrument)
#' b <- solve_b_pairs(projection_matrices(R, treatments))
#' target_populations(b)
#'
#' @export
solve_b_pairs <- function(projections, tolerance = 3) {
  treatments <- names(projections)
  pairs <- utils::combn(seq_along(treatments), 2L, simplify = FALSE)
  out <- lapply(pairs, function(ij) {
    solve_b_stack(list(projections[[ij[1]]]$K, projections[[ij[2]]]$K),
                  tolerance)
  })
  names(out) <- vapply(pairs, function(ij) {
    paste(treatments[ij[1]], treatments[ij[2]], sep = "_")
  }, character(1L))
  out
}

#' Describe the identified target populations
#'
#' Translates the binary solutions `b` into the adherence sets they select, so
#' that the target population of each identified effect can be read off
#' directly.
#'
#' @param b Output of [solve_b_pairs()] or [solve_b_treatments()].
#'
#' @return A list with the same names as `b`. Each element is a list of
#'   character vectors, one per solution, naming the adherence sets `a` for
#'   which `b[a] = 1`.
#'
#' @export
target_populations <- function(b) {
  lapply(b, function(bt) {
    apply(bt, 1L, function(row) names(which(row == 1)), simplify = FALSE)
  })
}
