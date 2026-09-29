# Build the pseudo-population behind an identified effect

Attaches a weight to each observation receiving one of the two
treatments in a contrast, such that a weighted comparison of outcomes
reproduces the estimate from
[`estimate_late()`](https://kazemia.github.io/catIV/reference/estimate_late.md).
The weighted sample is the target trial the instrumental variable
analysis emulates, so its weighted baseline characteristics describe the
population the effect applies to.

## Usage

``` r
pseudo_population(
  contrast,
  data,
  projections,
  b,
  p_sigma,
  instrument,
  treatment,
  solution = 1L
)
```

## Arguments

- contrast:

  A pair of treatments as a string, `"t1_t2"`.

- data:

  A data frame.

- projections:

  Output of
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

- b:

  Output of
  [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md).

- p_sigma:

  Output of
  [`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md).

- instrument, treatment:

  Column names, as strings.

- solution:

  Which identified target population to build the pseudo-population for,
  given as a solution name such as `"A4+A7"`, or as a positive integer
  indexing the rows of `b[[contrast]]`. Defaults to the first.

## Value

The subset of `data` receiving one of the two treatments, with an added
numeric column `w` holding the weights.

## Details

`p_sigma` must be the unadjusted estimate, computed from the same
[`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md)
output as the data being weighted. The weights are built from raw cell
proportions, so passing a covariate-standardised `P_Sigma` breaks the
correspondence: the weighted comparison no longer reproduces the
estimate exactly and the weights no longer average to one.

Given the unadjusted `p_sigma`, the weights average to exactly one
within each treatment arm. Both `weighted.mean(x, w)` and `mean(x * w)`
therefore give the same weighted average of a covariate within an arm,
which is what makes the weighted baseline table meaningful.

Weights can be negative. `b B_t^+` is a contrast across instrument
values and need not be non-negative, so some observations enter with a
negative sign. This is expected, not a numerical failure.

Called `PseudoPopulator()` in the paper code, where the target
population was selected by the integer `Pi_index`.
