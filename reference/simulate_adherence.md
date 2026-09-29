# Simulate from the adherence set model

Generates data the way the identification results assume it arises: a
confounder determines the decision team's adherence set, the instrument
is drawn independently, and categorical monotonicity then picks the
treatment deterministically. Effects estimated from this process are
exactly the quantities the package identifies, so it is the right
generator for checking that an estimator recovers a known truth.

## Usage

``` r
simulate_adherence(
  n,
  instrument,
  sets,
  treatments,
  v_levels,
  v_probs,
  v_to_set,
  v_on_y,
  t_on_y,
  choice = choose_treatment
)
```

## Arguments

- n:

  Number of observations.

- instrument:

  A named list of instrument values, as for
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).
  Drawn uniformly.

- sets:

  A named list of adherence sets, as for
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).

- treatments:

  Character vector of the treatment alternatives.

- v_levels:

  Values the confounder can take.

- v_probs:

  Probability of each confounder level; must sum to one.

- v_to_set:

  A list with one element per confounder level, named by that level,
  each a vector of probabilities over `sets` summing to one.

- v_on_y:

  Effect of each confounder level on the outcome probability, in the
  order of `v_levels`.

- t_on_y:

  Effect of each treatment on the outcome probability, in the order of
  `treatments`.

- choice:

  The choice function, by default
  [`choose_treatment()`](https://kazemia.github.io/catIV/reference/choose_treatment.md).

## Value

A data frame with columns `Z` (the name of the instrument value, or its
index if `instrument` is unnamed), `T`, `Y`, and `A` naming the
adherence set each observation was drawn with. `A` is unobservable in
practice and is returned so that simulations can check against it.

## Details

The outcome is Bernoulli with probability `t_on_y[T] + v_on_y[V]`, so
the caller must choose effects keeping every total in `[0, 1]`.

Called `CatSimulator()` in paper 1's code, which did not return the
adherence set.
