# Conditional local average treatment responses

The categorical instrumental variable estimator of paper 3: for each
treatment and each solution `b`, estimates
`E(Y^{T=t} | b[A] = 1, X = x)`.

## Usage

``` r
conditional_latr(
  x,
  numerator_model,
  denominator_model,
  projections,
  b,
  instrument_levels,
  treatments,
  instrument = "Z",
  treatment = "T",
  covariate = "X"
)
```

## Arguments

- x:

  A single value of the effect-modifying covariate.

- numerator_model:

  A fitted model of the outcome on the instrument, the treatment and the
  covariate, for example `lm(y ~ Z * trt * X, data = d)`.

- denominator_model:

  A fitted model of the treatment on the instrument and the covariate,
  as for
  [`conditional_p_z()`](https://kazemia.github.io/catIV/reference/conditional_p_z.md).

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

- treatment:

  Name of the treatment column as the numerator model was fitted with
  it.

## Value

A named list over treatments of named numeric vectors over the solutions
in `b`.

## Details

The numerator `Q_Z(t, X)` is formed by multiplying the modelled outcome
means `E(Y | Z, T, X)` by the modelled treatment probabilities
`P(T | Z, X)`, then projecting through `b B_t^+`, exactly as in the
unconditional case. Neither model is a local estimator, so the variance
expression in the paper is an approximation here.

Called `CIV_estimator()` in the paper 3 code.
