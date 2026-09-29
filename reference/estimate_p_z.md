# Estimate the conditional treatment probabilities P_Z

`P_Z(t)[i] = P(T = t | Z = z_i)`, the probability of receiving each
treatment under each value of the instrument.

## Usage

``` r
estimate_p_z(
  data,
  instrument,
  treatment,
  covariates = NULL,
  parametric = FALSE,
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

- covariates:

  Character vector of columns to adjust for. Only available when
  `parametric = TRUE`; supplying them otherwise is an error rather than
  being silently ignored.

- parametric:

  Estimate by multinomial regression rather than by cell proportions.
  Required for covariate adjustment.

- instrument_levels, treatments:

  Optional character vectors fixing the row and column order. Default to
  the sorted unique values found in `data`. Pass them explicitly to
  guarantee the rows line up with the rows of
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).

## Value

A numeric matrix with one row per instrument value and one column per
treatment. Rows sum to one.

## Details

The parametric estimator fits a multinomial regression of treatment on
the instrument and the covariates, then averages the predicted
probabilities over the observed covariate distribution separately at
each value of the instrument. That standardisation is what makes the
adjusted estimator a valid input to the identification results.

Called `MakeP_Z()` in the paper code.

## See also

[`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md),
[`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md)

## Examples

``` r
d <- data.frame(
  z = rep(1:3, each = 40),
  trt = rep(c("a", "b", "a", "b", "b", "a"), each = 20)
)
estimate_p_z(d, "z", "trt")
#>     a   b
#> 1 0.5 0.5
#> 2 0.5 0.5
#> 3 0.5 0.5
```
