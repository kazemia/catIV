# A naive instrumental variable estimate ignoring the other treatments

Collapses the ordinal instrument to a binary one, indicating whether the
first treatment is more encouraged than the second, restricts the data
to observations receiving one of the two, and applies the usual Wald
ratio.

## Usage

``` r
naive_iv(contrast, data, instrument_values, instrument, treatment, outcome)
```

## Arguments

- contrast:

  A pair of treatments as a string, `"t1_t2"`.

- data:

  A data frame.

- instrument_values:

  The named list of instrument values, as passed to
  [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md).

- instrument, treatment, outcome:

  Column names, as strings.

## Value

A single number: the Wald estimate.

## Details

This is the comparison estimator in paper 1, showing what happens when
the remaining treatment alternatives are ignored rather than modelled.
It is generally biased under categorical monotonicity, because switching
between the two treatments of interest is not the only response to the
instrument.

Called `NaiveIV()` in the paper code.
