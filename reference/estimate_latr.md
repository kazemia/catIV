# Identify local average treatment responses

The single-arm counterpart of
[`estimate_late()`](https://kazemia.github.io/catIV/reference/estimate_late.md):
for each treatment and each solution `b`, estimates
`E(Y^{T=t} | b[A] = 1)`.

## Usage

``` r
estimate_latr(q_z, projections, b, p_sigma)
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
  [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md),
  indexed by single treatments rather than by contrasts.

- p_sigma:

  A named list over treatments of named numeric vectors, as returned by
  [`estimate_p_sigma_treatments()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma_treatments.md).

## Value

A named list over treatments of named numeric vectors over the solutions
in `b`.

## Details

Called `LATOIdentifier()` in the paper code.
