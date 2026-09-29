# Identify the target populations of local average treatment effects

Finds every binary vector `b` that solves
`b = b B_t^+ B_t = b B_t'^+ B_t'` simultaneously, for each pair of
treatment alternatives. Each solution defines a target population for
which the contrast `E(Y^{T=t} - Y^{T=t'} | b[A] = 1)` is identified.

## Usage

``` r
solve_b_pairs(projections, tolerance = 3)
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

A named list with one matrix per pair of treatments, named `"t1_t2"`.
Each row is a solution `b`, each column an adherence set. A pair with no
solution yields a zero-row matrix, meaning no effect is identified for
that contrast.

## Examples

``` r
instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"),
                   z3 = c("b", "c", "a"), z4 = c("c", "b", "a"))
treatments <- c("a", "b", "c")
R <- response_matrix(adherence_sets(treatments), instrument)
b <- solve_b_pairs(projection_matrices(R, treatments))
target_populations(b)
#> $a_b
#> $a_b$`A4+A7`
#> [1] "A4" "A7"
#> 
#> 
#> $a_c
#> $a_c$A5
#> [1] "A5"
#> 
#> 
#> $b_c
#> $b_c$`A6+A7`
#> [1] "A6" "A7"
#> 
#> 
```
