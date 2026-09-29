#' catIV: Instrumental Variable Analysis with a Categorical Treatment
#'
#' Identification and estimation of causal effects when an ordinal instrument
#' affects the choice of a categorical treatment, following the adherence set
#' framework and the categorical monotonicity assumption.
#'
#' @section Where to start:
#' A typical analysis runs through four stages:
#'
#' 1. Describe the design. [adherence_sets()] enumerates the possible
#'    adherence sets, or supply a shorter list of clinically realistic ones.
#'    [response_matrix()] records which treatment each set receives under each
#'    value of the instrument.
#' 2. Work out what is identified. [projection_matrices()] builds `B_t` and
#'    `K_t`; [solve_b_pairs()] finds the target populations for which a
#'    contrast is identified, and [target_populations()] says which adherence
#'    sets each one contains.
#' 3. Estimate the observable quantities from data with [estimate_p_z()] and
#'    [estimate_q_z()], with or without covariate adjustment.
#' 4. Combine them. [estimate_p_sigma()] gives the size of each target
#'    population and [estimate_late()] the local average treatment effects.
#'
#' @section Function names used in the papers:
#' The published code for the three papers used different names. The mapping
#' is:
#'
#' \tabular{lll}{
#'   **Paper code**        \tab **catIV**                       \tab **Notes** \cr
#'   `GenerateA()`         \tab [adherence_sets()]              \tab \cr
#'   `T_decider()`         \tab [choose_treatment()]            \tab \cr
#'   `MakeR()`             \tab [response_matrix()]             \tab \cr
#'   `MakeKB()`            \tab [projection_matrices()]         \tab list elements `B_t`, `B_t_i` are now `B`, `B_plus` \cr
#'   `KbSolver()`          \tab [solve_b_pairs()]               \tab papers 1 and 2, indexed by treatment pairs \cr
#'   `KbSolver()`          \tab [solve_b_treatments()]          \tab paper 3, indexed by single treatments \cr
#'   `PiIdentifier()`      \tab [target_populations()]          \tab \cr
#'   `MakeP_Z()`           \tab [estimate_p_z()]                \tab \cr
#'   `MakeQ_Z()`           \tab [estimate_q_z()]                \tab \cr
#'   `MakeV_Z()`           \tab [estimate_v_z()]                \tab \cr
#'   `P_SigmaIdentifier()` \tab [estimate_p_sigma()]            \tab `E1`, `E2`, `WA` are now `arm1`, `arm2`, `average` \cr
#'   `P_SigmaIdentifier()` \tab [estimate_p_sigma_treatments()] \tab the paper 3 variant \cr
#'   `LATEIdentifier()`    \tab [estimate_late()]               \tab `RR` and `AverageProb` are now `scale` and `denominator` \cr
#'   `LATOIdentifier()`    \tab [estimate_latr()]               \tab
#' }
#'
#' @section Conventions:
#'
#' * `P_Z` and `Q_Z` are matrices with the instrument values as row names, so
#'   nothing depends on column position.
#' * The solvers always return a matrix. A contrast with no identifiable
#'   effect gives a matrix with zero rows.
#' * Solutions are named by the adherence sets they select, for example
#'   `"A4+A7"`, so results are addressed by meaning rather than by position.
#' * Every instrument value must be observed in the data, so that the rows of
#'   `P_Z` line up with the rows of the response matrix.
#'
#' @keywords internal
#' @importFrom foreach %dopar%
"_PACKAGE"
