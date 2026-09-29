# Enumerate all possible adherence sets

An adherence set is the set of treatment alternatives that a decision
team would consider for a given patient, that is the alternatives with
non-zero probability of being chosen once the effect of the instrument
is marginalised out. With `NT` treatment alternatives there are
`2^NT - 1` non-empty adherence sets.

## Usage

``` r
adherence_sets(treatments)
```

## Arguments

- treatments:

  Character vector of the treatment alternatives, `Supp(T)`.

## Value

A named list of character vectors, one per adherence set, named `A1`,
`A2`, ... in order of increasing set size.

## Details

In applications the full enumeration is often implausible, and a smaller
list of clinically realistic adherence sets is supplied to
[`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md)
instead. Restricting the list makes more effects identifiable, at the
cost of an additional substantive assumption.

## Examples

``` r
adherence_sets(c("a", "b", "c"))
#> $A1
#> [1] "a"
#> 
#> $A2
#> [1] "b"
#> 
#> $A3
#> [1] "c"
#> 
#> $A4
#> [1] "a" "b"
#> 
#> $A5
#> [1] "a" "c"
#> 
#> $A6
#> [1] "b" "c"
#> 
#> $A7
#> [1] "a" "b" "c"
#> 
```
