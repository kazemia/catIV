# Identify the target populations of local average treatment responses

Finds every binary vector `b` solving `b = b B_t^+ B_t` for each
treatment alternative separately. Each solution defines a target
population, namely the decision teams whose adherence set `a` satisfies
`b[a] = 1`, for which the local average treatment response
`E(Y^{T=t} | b[A] = 1)` is identified.

## Usage

``` r
solve_b_treatments(projections, tolerance = 3)
```

## Arguments

- projections:

  Output of
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

- tolerance:

  Number of decimal places of accuracy. Must be *smaller* than the
  `tolerance` used in
  [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md).

## Value

A named list with one matrix per treatment. Each row is a solution `b`,
each column an adherence set.

## Details

This is the single-arm quantity used by the conditional estimators. For
contrasts between two treatments use
[`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md).

## Examples

``` r
instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"),
                   z3 = c("b", "c", "a"), z4 = c("c", "b", "a"))
treatments <- c("a", "b", "c")
R <- response_matrix(adherence_sets(treatments), instrument)
solve_b_treatments(projection_matrices(R, treatments))
#> $a
#>             A1 A2 A3 A4 A5 A6 A7
#> A5           0  0  0  0  1  0  0
#> A1           1  0  0  0  0  0  0
#> A4+A7        0  0  0  1  0  0  1
#> A1+A5        1  0  0  0  1  0  0
#> A4+A5+A7     0  0  0  1  1  0  1
#> A1+A4+A7     1  0  0  1  0  0  1
#> A1+A4+A5+A7  1  0  0  1  1  0  1
#> 
#> $b
#>             A1 A2 A3 A4 A5 A6 A7
#> A6+A7        0  0  0  0  0  1  1
#> A4+A7        0  0  0  1  0  0  1
#> A2+A6        0  1  0  0  0  1  0
#> A2+A4        0  1  0  1  0  0  0
#> A2+A4+A6+A7  0  1  0  1  0  1  1
#> 
#> $c
#>             A1 A2 A3 A4 A5 A6 A7
#> A5           0  0  0  0  1  0  0
#> A3           0  0  1  0  0  0  0
#> A6+A7        0  0  0  0  0  1  1
#> A3+A5        0  0  1  0  1  0  0
#> A5+A6+A7     0  0  0  0  1  1  1
#> A3+A6+A7     0  0  1  0  0  1  1
#> A3+A5+A6+A7  0  0  1  0  1  1  1
#> 
```
