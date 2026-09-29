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

## Function names used in the papers

The published code for the three papers used different names. The
mapping is:

|  |  |  |
|----|----|----|
| **Paper code** | **catIV** | **Notes** |
| `GenerateA()` | [`adherence_sets()`](https://kazemia.github.io/catIV/reference/adherence_sets.md) |  |
| `T_decider()` | [`choose_treatment()`](https://kazemia.github.io/catIV/reference/choose_treatment.md) |  |
| `MakeR()` | [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md) |  |
| `MakeKB()` | [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md) | list elements `B_t`, `B_t_i` are now `B`, `B_plus` |
| `KbSolver()` | [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md) | papers 1 and 2, indexed by treatment pairs |
| `KbSolver()` | [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md) | paper 3, indexed by single treatments |
| `PiIdentifier()` | [`target_populations()`](https://kazemia.github.io/catIV/reference/target_populations.md) |  |
| `MakeP_Z()` | [`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md) |  |
| `MakeQ_Z()` | [`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md) |  |
| `MakeV_Z()` | [`estimate_v_z()`](https://kazemia.github.io/catIV/reference/estimate_v_z.md) |  |
| `P_SigmaIdentifier()` | [`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md) | `E1`, `E2`, `WA` are now `arm1`, `arm2`, `average` |
| `P_SigmaIdentifier()` | [`estimate_p_sigma_treatments()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma_treatments.md) | the paper 3 variant |
| `LATEIdentifier()` | [`estimate_late()`](https://kazemia.github.io/catIV/reference/estimate_late.md) | `RR` and `AverageProb` are now `scale` and `denominator` |
| `LATOIdentifier()` | [`estimate_latr()`](https://kazemia.github.io/catIV/reference/estimate_latr.md) |  |

## Differences from the paper code

The results are numerically identical to the published code, but a few
conventions changed:

- `P_Z` and `Q_Z` are returned as matrices with the instrument values as
  row names, rather than as data frames whose first column holds the
  instrument. Nothing now depends on column position.

- The solvers always return a matrix. A contrast with no identifiable
  effect gives a matrix with zero rows rather than a bare vector.

- Solutions are named by the adherence sets they select, for example
  `"A4+A7"`, so results are addressed by meaning rather than by
  position. Solution order therefore differs from the paper code, and
  index-based references such as `Pi_index = 1` do not carry over.

- Every instrument value must be observed in the data. The paper code
  silently produced a shorter `P_Z` in that situation.

## See also

Useful links:

- <https://github.com/kazemia/catIV>

- <https://kazemia.github.io/catIV/>

- Report bugs at <https://github.com/kazemia/catIV/issues>

## Author

**Maintainer**: Amir Aamodt Kazemi <amirhosk@uio.no>
