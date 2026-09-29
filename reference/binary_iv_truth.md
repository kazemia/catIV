# True principal stratum probabilities for a binary instrument

Computes the exact conditional probabilities of each principal stratum,
and of receiving treatment, for the binary data-generating process used
in the simulations of paper 3. Because these are the values the
estimators are trying to recover, they give a simulation something to be
checked against.

## Usage

``` r
binary_iv_truth(
  x,
  z_on_t,
  x_on_t,
  u_on_t,
  v_on_t,
  u_prob,
  v_prob,
  z_prob = 0.5
)
```

## Arguments

- x:

  Values of the effect-modifying covariate.

- z_on_t, x_on_t, u_on_t, v_on_t:

  Log-odds effects of the instrument, the covariate, the unobserved
  confounder and the observed confounder on the treatment.

- u_prob, v_prob, z_prob:

  Prevalences of the unobserved confounder, the observed confounder and
  the instrument.

## Value

A data frame with one row per value of `x`:

- `complier`:

  `P(T = 1 | Z = 1, X) - P(T = 1 | Z = 0, X)`.

- `always_taker`:

  `P(T = 1 | Z = 0, X)`.

- `never_taker`:

  `P(T = 0 | Z = 1, X)`.

- `treated`:

  `P(T = 1 | X)`, marginal over the instrument.

The three stratum probabilities sum to one.

## Details

The process is
`P(T = 1 | Z, U, V, X) = plogis(z_on_t * Z + u_on_t * U + v_on_t * V + x_on_t * X)`
with `U`, `V` and `Z` independent Bernoulli variables. `U` is unobserved
confounding, `V` is observed confounding.

This replaces `C_probability()`, `AT_probability()`, `NT_probability()`
and `T_probability()` from the paper 3 code. Those were written as four
separate expressions, and `T_probability()` took no `VT` or `VP`
arguments, so it marginalised over the unobserved confounder and the
instrument but not over the observed confounder. Setting `v_on_t = 0`
here reproduces that behaviour; any other value corrects it. Everything
is marginalised over the same 2x2 grid of `U` and `V`, so the four
quantities cannot drift apart.

## Examples

``` r
binary_iv_truth(c(0, 1, 2), z_on_t = 1.5, x_on_t = 0.4, u_on_t = 1,
                v_on_t = 0.8, u_prob = 0.5, v_prob = 0.5, z_prob = 0.5)
#>   x  complier always_taker never_taker   treated
#> 1 0 0.2089600    0.6947955  0.09624446 0.7992755
#> 2 1 0.1652479    0.7674115  0.06734060 0.8500354
#> 3 2 0.1257232    0.8277508  0.04652598 0.8906124
```
