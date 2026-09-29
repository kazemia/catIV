# Bootstrap confidence intervals for local average treatment effects

Resamples the data with replacement, re-estimates `P_Z`, `Q_Z`,
`P_Sigma` and the effects on each resample, and takes empirical
quantiles of the resulting distribution.

## Usage

``` r
bootstrap_late(
  n,
  data,
  instrument,
  treatment,
  outcome,
  projections,
  b,
  alpha = 0.05,
  n_cores = 14,
  cap = FALSE,
  covariates = NULL,
  parametric = FALSE,
  family = NULL,
  seed = 1234,
  instrument_levels = NULL,
  treatments = NULL
)
```

## Arguments

- n:

  Number of bootstrap replicates. With multiply imputed data this is the
  number drawn from *each* imputation.

- data:

  A data frame, or a list of data frames holding multiple imputations of
  the same study. With a list, `n` replicates are drawn from each
  imputation and all replicates are pooled before taking quantiles.

- instrument, treatment, outcome:

  Column names, as strings.

- projections:

  Output of
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

- b:

  Output of
  [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md).

- alpha:

  Two-sided error rate; the interval runs from the `alpha / 2` to the
  `1 - alpha / 2` quantile.

- n_cores:

  Number of worker processes. **This affects the numbers**: the
  replicates are split into one chunk per worker, so the random number
  streams, and hence the draws, depend on it. Fix it to reproduce a
  previous run. The published analyses used 14.

- cap:

  If `TRUE`, discard replicates whose estimated effect lies outside the
  observed range of the outcome before taking quantiles. Ratio
  estimators with a small denominator occasionally produce extreme
  values, and those replicates would otherwise dominate the quantiles.

- covariates, parametric, family:

  Passed to
  [`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md)
  and
  [`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md)
  to bootstrap the covariate-adjusted estimator.

- seed:

  Seed for the parallel random number streams.

- instrument_levels, treatments:

  Passed to the estimators to fix the row and column order. Default to
  the dimnames of the projection matrices.

## Value

A list with three elements:

- `alpha`:

  the error rate used.

- `replicates`:

  a data frame with one row per bootstrap replicate and one column per
  identified effect.

- `ci`:

  a data frame of `lower` and `upper` bounds, plus the number of
  replicates each interval is based on.

## Details

Replicates in which an instrument value drops out of the resample fail
to estimate, and are recorded as `NA` rather than aborting the run. The
`n_kept` column of `ci` reports how many replicates each interval
actually used, so that a badly behaved contrast is visible rather than
silent.

Called `BSCICalculator()` in the paper code, where the
multiple-imputation case was selected with `Data.complete = FALSE`.

## Examples

``` r
# \donttest{
treatments <- c("a", "b")
instrument <- list(`1` = c("a", "b"), `2` = c("b", "a"))
R <- response_matrix(adherence_sets(treatments), instrument)
KB <- projection_matrices(R, treatments)
b <- solve_b_pairs(KB)
d <- data.frame(
  z = rep(1:2, each = 200),
  trt = c(sample(c("a", "b"), 200, TRUE, c(0.7, 0.3)),
          sample(c("a", "b"), 200, TRUE, c(0.3, 0.7))),
  y = rbinom(400, 1, 0.5)
)
bootstrap_late(20, d, "z", "trt", "y", KB, b, n_cores = 1)$ci
#>            lower     upper n_kept
#> a_b.A3 -0.420505 0.1950215     20
# }
```
