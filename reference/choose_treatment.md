# The choice function implied by categorical monotonicity

Categorical monotonicity states that the treatment received is a
deterministic function of the instrument and the adherence set: the
decision team picks the alternative in the adherence set that the
instrument encourages most.

## Usage

``` r
choose_treatment(z, a)
```

## Arguments

- z:

  Character vector of all treatment alternatives, sorted in *decreasing*
  order of encouragement by the instrument. In the pricing application
  this is the alternatives ordered from cheapest to most expensive.

- a:

  Character vector of the alternatives the decision team adheres to,
  i.e. one adherence set.

## Value

A length-one character vector: the element of `a` appearing earliest in
`z`.

## Examples

``` r
choose_treatment(c("a", "b", "c"), c("b", "c"))
#> [1] "b"
```
