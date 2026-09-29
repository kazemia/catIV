# Heterogeneous instrumental variable estimator of conditional responses

The estimator of paper 3. A confounded outcome model is efficient but
biased by unobserved confounding; the categorical instrumental variable
estimator is unbiased but unstable wherever the target population is
small. This fits a weighted model for the difference between the two,
and subtracts the fitted bias from the confounded estimate, giving a
conditional average treatment response that keeps most of the efficiency
of the confounded model without its bias.

## Usage

``` r
hiv_cate(
  x_values,
  numerator_model,
  denominator_model,
  confounded_model,
  projections,
  b,
  density,
  instrument_levels,
  treatments,
  instrument = "Z",
  treatment = "T",
  covariate = "X",
  density_power = 1
)
```

## Arguments

- x_values:

  Numeric vector of covariate values to estimate at.

- numerator_model, denominator_model:

  Models of the outcome and of the treatment, as for
  [`conditional_latr()`](https://kazemia.github.io/catIV/reference/conditional_latr.md).

- confounded_model:

  An outcome model adjusting only for observed covariates, as for
  [`confounded_response()`](https://kazemia.github.io/catIV/reference/confounded_response.md).

- projections:

  Output of
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

- b:

  Output of
  [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md).

- density:

  Estimated density of the covariate at `x_values`, as a numeric vector
  of the same length. See details.

- instrument_levels, treatments:

  Character vectors of the instrument values and treatment alternatives.

- instrument, treatment, covariate:

  Column names as the models were fitted with them.

- density_power:

  Exponent on the density in the regression weights. Defaults to 1,
  matching the variance argument in the paper. Use 2 to reproduce the
  published analysis code. See the Weights section.

## Value

A data frame with one row per covariate value and treatment:

- `x`, `treatment`:

  where the estimate applies.

- `confounded`:

  the covariate-adjusted response.

- `bias`:

  the fitted bias.

- `response`:

  the debiased conditional average treatment response,
  `confounded - bias`.

The fitted bias model is attached as the `"bias_model"` attribute and
the assembled per-solution data as `"bias_data"`.

## Details

This estimator is assembled from the analysis script for paper 3 rather
than from its function file. Its components,
[`conditional_p_sigma()`](https://kazemia.github.io/catIV/reference/conditional_p_sigma.md)
and
[`conditional_latr()`](https://kazemia.github.io/catIV/reference/conditional_latr.md),
reproduce the originals exactly.

Inference is by bootstrapping the whole procedure; there is no closed
form.

## Weights

Each observation of the bias enters the weighted regression with weight
`density^density_power * p_sigma^2`, where `p_sigma` is the conditional
probability of belonging to the target population, that is the
denominator of the instrumental variable estimator.

The paper specifies the weight as the square of that denominator times
the density of the covariate, consistent with its variance expression,
in which the variance is proportional to `1 / (p_sigma^2 * density)`.
That gives `density_power = 1`, which is the default here.

The analysis code used for the published clinical results computed the
weight with the density **squared**. Set `density_power = 2` to
reproduce those numbers exactly. The default does not reproduce them,
because the squared density was not what the paper specifies.

## Contrasts

Responses are estimated one treatment at a time. Form a conditional
average treatment effect by differencing two of them, which
[`hiv_contrast()`](https://kazemia.github.io/catIV/reference/hiv_contrast.md)
does.
