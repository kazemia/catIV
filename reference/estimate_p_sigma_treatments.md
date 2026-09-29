# Probability of belonging to each single-treatment target population

The single-arm counterpart of
[`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md),
for use with
[`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md)
and
[`estimate_latr()`](https://kazemia.github.io/catIV/reference/estimate_latr.md).

## Usage

``` r
estimate_p_sigma_treatments(p_z, projections, b)
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
  [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md).

## Value

A named list over treatments of named numeric vectors over the solutions
in `b`.

## Details

This is the version of `P_SigmaIdentifier()` defined in `HTEfunctions.R`
for paper 3, which indexes `b` by single treatments.
