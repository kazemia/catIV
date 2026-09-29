# Adjustment, inference and target trials

``` r

library(catIV)
```

This vignette picks up where
[`vignette("identification")`](https://kazemia.github.io/catIV/articles/identification.md)
leaves off: you have an identified effect, and now you need a confidence
interval for it and a way to say who it applies to. Both matter more
here than in a randomised trial, because the target population is not
the enrolled sample but an unobservable subgroup defined by adherence
sets.

## Setup

``` r

treatments <- c("cheap", "mid", "dear")
instrument <- list(
  "2019" = c("cheap", "mid",   "dear"),
  "2020" = c("mid",   "cheap", "dear"),
  "2021" = c("mid",   "dear",  "cheap"),
  "2022" = c("dear",  "mid",   "cheap")
)

sets <- adherence_sets(treatments)
R    <- response_matrix(sets, instrument)
KB   <- projection_matrices(R, treatments)
b    <- solve_b_pairs(KB)
```

Simulated data with a measured baseline covariate `age` that predicts
both treatment and outcome:

``` r

set.seed(2)
d <- simulate_adherence(
  n = 6000, instrument = instrument, sets = sets, treatments = treatments,
  v_levels = c(0, 1), v_probs = c(0.5, 0.5),
  v_to_set = list("0" = c(3, 1, 1, 2, 4, 1, 3),
                  "1" = c(1, 3, 2, 1, 2, 3, 4)),
  v_on_y = c(0.10, 0.30), t_on_y = c(0.20, 0.35, 0.50)
)
d$age <- round(50 + 8 * (d$A %in% c("A3", "A6", "A7")) + rnorm(nrow(d), 0, 6))
```

## Covariate adjustment

The unadjusted estimators use raw cell proportions. The adjusted
versions fit a model and then **standardise** over the observed
covariate distribution at each instrument value, which is what keeps
them valid inputs to the identification algebra. A per-covariate-profile
prediction would not be.

``` r

p_z_adj <- estimate_p_z(d, "Z", "T", covariates = "age", parametric = TRUE,
                        instrument_levels = rownames(R),
                        treatments = treatments)
q_z_adj <- estimate_q_z(d, "Z", "T", "Y", covariates = "age",
                        parametric = TRUE, family = "binomial",
                        p_z = p_z_adj, instrument_levels = rownames(R),
                        treatments = treatments)

p_sigma_adj <- estimate_p_sigma(p_z_adj, KB, b)
lapply(estimate_late(q_z_adj, KB, b, p_sigma_adj), round, 3)
#> $cheap_mid
#>  A4+A7 
#> -0.132 
#> 
#> $cheap_dear
#>     A5 
#> -0.341 
#> 
#> $mid_dear
#>  A6+A7 
#> -0.121
```

Note that adjustment is refused unless you ask for it explicitly:

``` r

estimate_p_z(d, "Z", "T", covariates = "age")
#> Error:
#> ! Covariate adjustment requires `parametric = TRUE`.
```

Because the instrument is as good as randomly assigned here, adjustment
buys precision rather than correcting bias. It matters when the
instrument is only conditionally valid.

## Confidence intervals

There is no closed-form variance for these ratio estimators, so
inference is by bootstrap.
[`bootstrap_late()`](https://kazemia.github.io/catIV/reference/bootstrap_late.md)
resamples, re-estimates everything on each resample, and takes empirical
quantiles.

``` r

bs <- bootstrap_late(300, d, "Z", "T", "Y", KB, b, n_cores = 2)
round(bs$ci, 3)
#>                  lower  upper n_kept
#> cheap_mid.A4+A7 -0.226 -0.042    300
#> cheap_dear.A5   -0.441 -0.240    300
#> mid_dear.A6+A7  -0.207 -0.044    300
```

Three things worth knowing:

- **`n_cores` changes the numbers.** Replicates are split one chunk per
  worker, so the random number streams depend on how many workers there
  are. Fix `n_cores` and `seed` together to reproduce a run.
- **`n_kept` is not decoration.** A resample that loses an instrument
  value cannot be estimated and is recorded as `NA`. If `n_kept` is well
  below the number of replicates, the interval rests on less than it
  appears to.
- **`cap = TRUE`** discards replicates whose estimate falls outside the
  range of the outcome. When a target population is small, the
  denominator of the ratio occasionally lands near zero and throws out
  an absurd value that would otherwise drag the quantiles.

``` r

# With multiply imputed data, pass a list. Each imputation contributes n
# replicates and all of them are pooled before quantiles are taken.
bootstrap_late(300, list(imp1, imp2, imp3), "Z", "T", "Y", KB, b)
```

## Describing the target population

The awkward part of reporting these estimates is saying who they apply
to. “Decision teams whose adherence set contains both alternatives” is
correct but not a sentence a clinician can act on. The target trial
framing helps: if the adherence set *were* observable, you could run a
trial that enrolled the teams in the target population and randomised
them among the alternatives they were willing to accept. The IV analysis
estimates the result of that trial.

[`pseudo_population()`](https://kazemia.github.io/catIV/reference/pseudo_population.md)
builds the weighted dataset corresponding to that trial. Give it the
**unadjusted** `p_sigma`: the weights are built from the raw data, so
mixing in a standardised `P_Sigma` breaks the correspondence.

``` r

p_z     <- estimate_p_z(d, "Z", "T", instrument_levels = rownames(R),
                        treatments = treatments)
q_z     <- estimate_q_z(d, "Z", "T", "Y", instrument_levels = rownames(R),
                        treatments = treatments)
p_sigma <- estimate_p_sigma(p_z, KB, b)

pseudo <- pseudo_population("cheap_mid", d, KB, b, p_sigma, "Z", "T")
nrow(pseudo)
#> [1] 4331
range(pseudo$w)
#> [1] -5.021207  4.083763
```

Weights can be negative. `b B_t^+` is a contrast across instrument
values, not a probability, so some observations enter with a negative
sign. That is expected rather than a numerical failure.

Weighting makes an ordinary comparison of outcomes reproduce the IV
estimate:

``` r

arm1 <- pseudo$T == "cheap"
weighted.mean(pseudo$Y[arm1], pseudo$w[arm1]) -
  weighted.mean(pseudo$Y[!arm1], pseudo$w[!arm1])
#> [1] -0.133972

estimate_late(q_z, KB, b, p_sigma)$cheap_mid
#>     A4+A7 
#> -0.133972
```

The point of the pseudo-population is not the effect estimate, though —
it is the baseline table. Weighted covariate means describe the
population the effect applies to, in the same way a trial’s Table 1
describes who was enrolled:

``` r

tapply(seq_len(nrow(pseudo)), pseudo$T, function(i) {
  round(weighted.mean(pseudo$age[i], pseudo$w[i]), 1)
})
#> cheap   mid 
#>  55.1  55.7

# Unweighted, for comparison
round(tapply(d$age, d$T, mean)[c("cheap", "mid")], 1)
#> cheap   mid 
#>  51.1  54.2
```

The weights average to exactly one within each treatment arm, so
`mean(x * w)` within an arm equals the weighted mean of `x` and either
form gives the same table.

``` r

tapply(pseudo$w, pseudo$T, mean)
#> cheap   mid 
#>     1     1
```

## Choosing what to report

With several identified contrasts, not all deserve reporting. A contrast
whose target population is tiny will have a wide interval and describe
almost nobody.
[`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md)
gives the sizes, and screening on them is how paper 2 chose which
comparisons to present:

``` r

sizes <- vapply(p_sigma_adj$average, max, numeric(1))
round(sizes, 3)
#>  cheap_mid cheap_dear   mid_dear 
#>      0.316      0.206      0.352

names(sizes)[sizes > 0.1]
#> [1] "cheap_mid"  "cheap_dear" "mid_dear"
```
