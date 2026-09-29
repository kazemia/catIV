# Average a conditional effect over the covariate distribution

Collapses a conditional average treatment effect to a single number by
averaging over the covariate, weighted by its estimated density.

## Usage

``` r
average_contrast(contrast, density)
```

## Arguments

- contrast:

  Output of
  [`hiv_contrast()`](https://kazemia.github.io/catIV/reference/hiv_contrast.md).

- density:

  Estimated density of the covariate at `contrast$x`.

## Value

A single number.
