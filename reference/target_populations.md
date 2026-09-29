# Describe the identified target populations

Translates the binary solutions `b` into the adherence sets they select,
so that the target population of each identified effect can be read off
directly.

## Usage

``` r
target_populations(b)
```

## Arguments

- b:

  Output of
  [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md)
  or
  [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md).

## Value

A list with the same names as `b`. Each element is a list of character
vectors, one per solution, naming the adherence sets `a` for which
`b[a] = 1`.
