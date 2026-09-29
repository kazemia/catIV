# Estimate the outcome-probability products Q_Z

`Q_Z(t)[i] = E(Y | Z = z_i, T = t) * P(T = t | Z = z_i)`. Ratios of
`Q_Z` to `P_Z`, both projected through `b B_t^+`, give the local average
treatment responses.

## Usage

``` r
estimate_q_z(
  data,
  instrument,
  treatment,
  outcome,
  covariates = NULL,
  parametric = FALSE,
  family = NULL,
  p_z = NULL,
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

- covariates:

  Character vector of columns to adjust for. Only available when
  `parametric = TRUE`; supplying them otherwise is an error rather than
  being silently ignored.

- parametric:

  Estimate by multinomial regression rather than by cell proportions.
  Required for covariate adjustment.

- family:

  Family passed to [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
  for the parametric outcome model, for example `"binomial"` for a
  binary outcome.

- p_z:

  A matrix from
  [`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md).
  Required when `parametric = TRUE`, because the adjusted outcome means
  are multiplied by the adjusted treatment probabilities.

- instrument_levels, treatments:

  Optional character vectors fixing the row and column order. Default to
  the sorted unique values found in `data`. Pass them explicitly to
  guarantee the rows line up with the rows of
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).

## Value

A numeric matrix with one row per instrument value and one column per
treatment.

## Details

Without adjustment this reduces to the sum of the outcome within each
instrument-by-treatment cell, divided by the number of observations at
that value of the instrument. Cells with no observations contribute
zero, which is right here because `P(T = t | Z = z)` is zero for them
too.

Called `MakeQ_Z()` in the paper code.
