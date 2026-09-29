# Identify local average treatment effects

For each contrast and each solution `b`, estimates
`E(Y^{T=t} - Y^{T=t'} | b[A] = 1)`, the effect among decision teams
whose adherence set satisfies `b[A] = 1`.

## Usage

``` r
estimate_late(
  q_z,
  projections,
  b,
  p_sigma,
  scale = c("difference", "ratio"),
  denominator = c("arm", "average"),
  arms = FALSE
)
```

## Arguments

- q_z:

  A matrix from
  [`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md).

- projections:

  Output of
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

- b:

  Output of
  [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md).

- p_sigma:

  Output of
  [`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md).

- scale:

  `"difference"` for a difference in means, `"ratio"` for a risk ratio.

- denominator:

  Which target-population probability to divide the two arms by: `"arm"`
  uses each arm's own estimate, `"average"` uses the weighted average
  for both.

- arms:

  If `TRUE`, also return the two local average treatment *responses*
  that the contrast is built from.

## Value

With `arms = FALSE` (the default), a named list over contrasts of named
numeric vectors over the solutions in `b`. With `arms = TRUE`, a list of
three such objects: `arm1`, `arm2` and `effect`.

## Details

Note that with `scale = "ratio"` and `denominator = "average"` the
shared denominator cancels, so the result is the ratio of the raw
projected quantities.

Called `LATEIdentifier()` in the paper code, where `scale` and
`denominator` were the logical arguments `RR` and `AverageProb`.
