# catIV: Instrumental Variable Analysis with a Categorical Treatment

Identification and estimation of causal effects when an ordinal
instrument affects the choice of a categorical treatment, following the
adherence set framework and the categorical monotonicity assumption.

## Where to start

A typical analysis runs through four stages:

1.  Describe the design.
    [`adherence_sets()`](https://kazemia.github.io/catIV/reference/adherence_sets.md)
    enumerates the possible adherence sets, or supply a shorter list of
    clinically realistic ones.
    [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md)
    records which treatment each set receives under each value of the
    instrument.

2.  Work out what is identified.
    [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md)
    builds `B_t` and `K_t`;
    [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md)
    finds the target populations for which a contrast is identified, and
    [`target_populations()`](https://kazemia.github.io/catIV/reference/target_populations.md)
    says which adherence sets each one contains.

3.  Estimate the observable quantities from data with
    [`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md)
    and
    [`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md),
    with or without covariate adjustment.

4.  Combine them.
    [`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md)
    gives the size of each target population and
    [`estimate_late()`](https://kazemia.github.io/catIV/reference/estimate_late.md)
    the local average treatment effects.

## Conventions

- `P_Z` and `Q_Z` are matrices with the instrument values as row names,
  so nothing depends on column position.

- The solvers always return a matrix. A contrast with no identifiable
  effect gives a matrix with zero rows.

- Solutions are named by the adherence sets they select, for example
  `"A4+A7"`, so results are addressed by meaning rather than by
  position.

- Every instrument value must be observed in the data, so that the rows
  of `P_Z` line up with the rows of the response matrix.

## See also

Useful links:

- <https://github.com/kazemia/catIV>

- <https://kazemia.github.io/catIV/>

- Report bugs at <https://github.com/kazemia/catIV/issues>

## Author

**Maintainer**: Amir Aamodt Kazemi <amirhosk@uio.no>

Authors:

- Amir Aamodt Kazemi <amirhosk@uio.no>
