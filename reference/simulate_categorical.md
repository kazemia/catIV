# Simulate a categorical treatment chosen by a multinomial model

The generator used for the simulation studies in papers 1 and 2. Unlike
[`simulate_adherence()`](https://kazemia.github.io/catIV/reference/simulate_adherence.md),
the treatment is drawn from a multinomial model given the instrument and
several binary confounders, rather than through a deterministic choice
function. Categorical monotonicity therefore holds only approximately,
which is the point: it tests the estimators under a process that does
not exactly satisfy their assumptions.

## Usage

``` r
simulate_categorical(
  n,
  instrument,
  treatments,
  v_on_t,
  v_probs,
  v_on_y,
  t_on_y,
  intercept = 0,
  observed = TRUE,
  z_on_t = c("exp", "linear")
)
```

## Arguments

- n:

  Number of observations.

- instrument:

  A named list of instrument values, as for
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).
  Drawn uniformly.

- treatments:

  Character vector of the treatment alternatives.

- v_on_t:

  A matrix with one row per treatment and one column per confounder,
  giving the effect of each confounder on each treatment.

- v_probs:

  Prevalence of each confounder, one value per confounder.

- v_on_y:

  Effect of each confounder on the outcome, one value per confounder.

- t_on_y:

  Effect of each treatment on the outcome, in the order of `treatments`.

- intercept:

  Added to every outcome probability.

- observed:

  Logical vector saying which confounders are returned. A confounder set
  to `FALSE` acts as unobserved confounding. Recycled to the number of
  confounders, so the default returns all of them.

- z_on_t:

  `"exp"` to exponentiate the instrument's encouragement rank before
  combining, `"linear"` to use it directly.

## Value

A data frame with columns `Z` (the name of the instrument value, or its
index if `instrument` is unnamed), `T`, `Y`, and one column `V1`, `V2`,
... per observed confounder.

## Details

Outcome probabilities are clipped to `[0, 1]` after adding the treatment
effect, the confounder effects and the intercept.

Called `CatSimulator2()` in paper 1's code and `CatSimulator()` in paper
2's; paper 2's version is the same generator with `intercept = 0`.
