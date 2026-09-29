# Probability of belonging to each identified target population

For every solution `b`, `b B_t^+ P_Z(t)` is the probability that a
decision team's adherence set satisfies `b[A] = 1`. Because a contrast
involves two treatments, the probability can be computed from either
arm; the two agree in the population but not in a finite sample.

## Usage

``` r
estimate_p_sigma(p_z, projections, b)
```

## Arguments

- p_z:

  A matrix from
  [`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md).

- projections:

  Output of
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

- b:

  Output of
  [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md).

## Value

A list of three elements, each a named list over contrasts of named
numeric vectors over the solutions in `b`:

- `arm1`:

  computed from the first treatment of the pair.

- `arm2`:

  computed from the second treatment of the pair.

- `average`:

  the weighted average of the two, weighting each arm by its marginal
  share `colSums(p_z)`.

## Details

`average` is the quantity used to screen out contrasts whose target
population is too small to estimate reliably. Called
`P_SigmaIdentifier()` in the paper code, where the three elements were
named `E1`, `E2` and `WA`.
