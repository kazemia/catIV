# Contrast two treatments from a heterogeneous IV fit

Contrast two treatments from a heterogeneous IV fit

## Usage

``` r
hiv_contrast(fit, t1, t2)
```

## Arguments

- fit:

  Output of
  [`hiv_cate()`](https://kazemia.github.io/catIV/reference/hiv_cate.md).

- t1, t2:

  The two treatments to contrast, as strings. The result is the response
  under `t1` minus the response under `t2`.

## Value

A data frame with columns `x`, `contrast` and `estimate`.
