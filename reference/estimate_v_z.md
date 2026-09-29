# Estimate the conditional outcome variances V_Z

`V_Z(t)[i] = Var(Y | Z = z_i, T = t) * P(T = t | Z = z_i)^2`, the
contribution of each instrument-by-treatment cell to the variance of
[`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md).

## Usage

``` r
estimate_v_z(
  data,
  instrument,
  treatment,
  outcome,
  instrument_levels = NULL,
  treatments = NULL
)
```

## Arguments

- data:

  A data frame with one row per observation.

- instrument:

  Name of the instrument column, as a string.

- treatment:

  Name of the treatment column, as a string.

- outcome:

  Name of the outcome column, as a string.

- instrument_levels, treatments:

  Optional character vectors fixing the row and column order. Default to
  the sorted unique values found in `data`. Pass them explicitly to
  guarantee the rows line up with the rows of
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).

## Value

A numeric matrix with one row per instrument value and one column per
treatment. Cells holding fewer than two observations contribute zero,
since no within-cell variance can be estimated from them.

## Details

Called `MakeV_Z()` in the paper code. There is no parametric version.
