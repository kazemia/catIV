# Conditional probability of belonging to each target population

The conditional counterpart of
[`estimate_p_sigma_treatments()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma_treatments.md):
for each treatment and each solution `b`, the probability
`P(b[A] = 1 | X = x)` of a decision team being in the identified target
population at covariate value `x`.

## Usage

``` r
conditional_p_sigma(
  x,
  model,
  projections,
  b,
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

- projections:

  Output of
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

- b:

  Output of
  [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md).

- instrument_levels:

  Character vector of the instrument values, in the order used to build
  the projection matrices.

- treatments:

  Character vector of the treatment alternatives.

- instrument, covariate:

  Names of the instrument and covariate columns as the model was fitted
  with them.

## Value

A named list over treatments of named numeric vectors over the solutions
in `b`.

## Details

Positivity requires this to be strictly positive for the covariate
values of interest. Where it approaches zero the instrumental variable
estimator has a vanishing denominator and becomes unstable, which is the
problem the estimator in
[`hiv_cate()`](https://kazemia.github.io/catIV/reference/hiv_cate.md) is
designed to work around.

Called `Pi_prob_estimator()` in the paper 3 code.
