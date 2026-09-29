# Build the response matrix

The response matrix records which treatment each adherence set receives
under each value of the instrument. It is the design matrix from which
the indicator matrices `B_t` are formed.

## Usage

``` r
response_matrix(sets, instrument, choice = choose_treatment)
```

## Arguments

- sets:

  A named list of adherence sets, typically from
  [`adherence_sets()`](https://kazemia.github.io/catIV/reference/adherence_sets.md)
  or a hand-specified list of clinically realistic sets.

- instrument:

  A named list of instrument values. Each element is a character vector
  of all treatment alternatives sorted in decreasing order of
  encouragement.

- choice:

  The choice function. Defaults to
  [`choose_treatment()`](https://kazemia.github.io/catIV/reference/choose_treatment.md),
  the function implied by categorical monotonicity.

## Value

A data frame with one row per instrument value and one column per
adherence set, holding the treatment received.

## Examples

``` r
instrument <- list(z1 = c("a", "b", "c"), z2 = c("b", "a", "c"))
response_matrix(adherence_sets(c("a", "b", "c")), instrument)
#>    A1 A2 A3 A4 A5 A6 A7
#> z1  a  b  c  a  a  b  a
#> z2  a  b  c  b  a  b  b
```
