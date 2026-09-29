# Conditional treatment responses from a confounded outcome model

Reads `E(Y | T = t, X = x)` off an outcome model that adjusts only for
observed covariates. These estimates are efficient but biased by any
unobserved confounding;
[`hiv_cate()`](https://kazemia.github.io/catIV/reference/hiv_cate.md)
combines them with the instrumental variable estimates to remove that
bias.

## Usage

``` r
confounded_response(x, model, treatments, treatment = "T", covariate = "X")
```

## Arguments

- x:

  Values of the effect-modifying covariate.

- model:

  A fitted outcome model, for example a `glm` or a
  [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

- treatments:

  Character vector of the treatment alternatives.

- treatment, covariate:

  Names of the treatment and covariate columns as the model was fitted
  with them.

## Value

A data frame with columns `x`, `treatment` and `response`.

## Details

Covariates other than `x` are held at typical values, numeric ones at
their mean and categorical ones at their mode, as
[`marginaleffects::datagrid()`](https://rdrr.io/pkg/marginaleffects/man/datagrid.html)
does. They are not averaged over their distribution.

Called `Confounded_estimator()` and `CatConfounded_estimator()` in the
paper code.
