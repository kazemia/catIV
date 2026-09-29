# Simulate a binary treatment with a binary instrument

The two-treatment data-generating process from paper 1: a single binary
confounder affects both the treatment and the outcome, and a binary
instrument shifts the treatment only.

## Usage

``` r
simulate_binary(n, z_on_t, v_on_t, v_on_y, t_on_y, logit = FALSE)
```

## Arguments

- n:

  Number of observations.

- z_on_t:

  Effect of the instrument on the treatment.

- v_on_t, v_on_y:

  Effects of the confounder on the treatment and on the outcome.

- t_on_y:

  Effect of the treatment on the outcome.

- logit:

  If `FALSE` (the default) the effects are probabilities and combine
  additively, so the caller must keep every total in `[0, 1]`. If `TRUE`
  they are log odds and pass through a logistic link, which is safe for
  any values.

## Value

A data frame with columns `Z`, `T` and `Y`, all binary.

## Details

Called `BinarySimulator()` in the paper code, where `logit` was the
argument `OR`.

## Examples

``` r
set.seed(1)
head(simulate_binary(10, z_on_t = 0.4, v_on_t = 0.3,
                     v_on_y = 0.25, t_on_y = 0.2))
#>   Z T Y
#> 1 0 0 0
#> 2 0 0 0
#> 3 1 0 1
#> 4 0 0 0
#> 5 1 1 0
#> 6 0 0 0
```
