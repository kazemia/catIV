# Indicator matrices and null-space projectors

For each treatment alternative `t` this builds the binary indicator
matrix `B_t`, where `B_t[i, n]` is one when adherence set `n` receives
treatment `t` under instrument value `i`; its Moore-Penrose
pseudo-inverse `B_t^+`; and the projector `K_t = I - B_t^+ B_t` onto the
null space of `B_t`.

## Usage

``` r
projection_matrices(response, treatments, tolerance = 4)
```

## Arguments

- response:

  A response matrix from
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).

- treatments:

  Character vector of the treatment alternatives.

- tolerance:

  Number of decimal places to which `K_t` is rounded. This removes
  floating-point noise from the pseudo-inverse. It must be *larger* than
  the `tolerance` later passed to the solvers, so that entries rounded
  to zero here are also treated as zero there.

## Value

A named list with one element per treatment, each a list of `K`, `B` and
`B_plus`.

## Details

The vectors `b` that identify a local average treatment response are the
binary solutions of `b = b B_t^+ B_t`, equivalently `b K_t = 0`. See
[`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md)
and
[`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md).

## Examples

``` r
instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"),
                   z3 = c("b", "c", "a"), z4 = c("c", "b", "a"))
treatments <- c("a", "b", "c")
R <- response_matrix(adherence_sets(treatments), instrument)
KB <- projection_matrices(R, treatments)
KB[["a"]]$B
#>    A1 A2 A3 A4 A5 A6 A7
#> z1  1  0  0  1  1  0  1
#> z2  1  0  0  0  1  0  0
#> z3  1  0  0  0  0  0  0
#> z4  1  0  0  0  0  0  0
```
