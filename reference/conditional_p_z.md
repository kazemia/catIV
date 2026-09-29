# Conditional treatment probabilities from a fitted model

Reads `P(T = t | Z = z, X = x)` off a model of treatment given the
instrument and an effect-modifying covariate. This is the conditional
counterpart of
[`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md),
and the denominator of the categorical instrumental variable estimator.

## Usage

``` r
conditional_p_z(
  x,
  model,
  instrument_levels,
  treatments,
  instrument = "Z",
  covariate = "X"
)
```

## Arguments

- x:

  A single value of the effect-modifying covariate.

- model:

  A fitted model of treatment on the instrument and the covariate, such
  as [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) or
  a random forest. Anything whose
  [`predict()`](https://rdrr.io/r/stats/predict.html) method returns
  class probabilities will work.

- instrument_levels:

  Character vector of the instrument values, in the order used to build
  the projection matrices.

- treatments:

  Character vector of the treatment alternatives.

- instrument, covariate:

  Names of the instrument and covariate columns as the model was fitted
  with them.

## Value

A numeric matrix with one row per instrument value and one column per
treatment, the same shape as
[`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md)
returns.

## Details

A two-class model whose
[`predict()`](https://rdrr.io/r/stats/predict.html) returns a single
vector of probabilities is expanded to two columns, taking the vector to
be the probability of the second treatment.
