# Identification with a categorical treatment

``` r

library(catIV)
```

## The problem

With a binary treatment, an instrument identifies the effect among
*compliers*: those who take the treatment when encouraged and not
otherwise. With three or more alternatives that description falls apart.
A decision team might be willing to prescribe A or B but never C.
Another might consider B or C but never A. An instrument that shifts the
ranking of all three moves these two groups in different directions, and
“complier” no longer names a single population.

`catIV` handles this by making the set of alternatives a decision team
would actually consider — its **adherence set** — the central object.

## A price-ranking instrument

Three treatments, and a procurement ranking that changes each year. Each
instrument value is the alternatives sorted from most to least
encouraged, which in this setting means cheapest first.

``` r

treatments <- c("cheap", "mid", "dear")

instrument <- list(
  "2019" = c("cheap", "mid",   "dear"),
  "2020" = c("mid",   "cheap", "dear"),
  "2021" = c("mid",   "dear",  "cheap"),
  "2022" = c("dear",  "mid",   "cheap")
)
```

The behavioural assumption, **categorical monotonicity**, is that a
decision team takes whichever alternative in its adherence set the
instrument encourages most:

``` r

choose_treatment(instrument[["2019"]], c("mid", "dear"))
#> [1] "mid"
choose_treatment(instrument[["2022"]], c("mid", "dear"))
#> [1] "dear"
```

A team open to `mid` or `dear` takes `mid` in 2019, when `mid` is
cheaper, and `dear` in 2022, when the ranking flips. A team open only to
`dear` takes `dear` in both years and contributes nothing to
identification.

## Enumerating adherence sets

With three alternatives there are `2^3 - 1 = 7` non-empty adherence
sets.

``` r

sets <- adherence_sets(treatments)
sets
#> $A1
#> [1] "cheap"
#> 
#> $A2
#> [1] "mid"
#> 
#> $A3
#> [1] "dear"
#> 
#> $A4
#> [1] "cheap" "mid"  
#> 
#> $A5
#> [1] "cheap" "dear" 
#> 
#> $A6
#> [1] "mid"  "dear"
#> 
#> $A7
#> [1] "cheap" "mid"   "dear"
```

The **response matrix** records what each set receives under each
instrument value. It is the whole behavioural model in one table.

``` r

R <- response_matrix(sets, instrument)
R
#>         A1  A2   A3    A4    A5   A6    A7
#> 2019 cheap mid dear cheap cheap  mid cheap
#> 2020 cheap mid dear   mid cheap  mid   mid
#> 2021 cheap mid dear   mid  dear  mid   mid
#> 2022 cheap mid dear   mid  dear dear  dear
```

Rows `A1` to `A3` are constant down every column: a team considering
only one alternative always takes it, whatever the price. Those teams
carry no information. The variation in the other columns is what
identifies effects.

## What is identifiable, and for whom?

From the response matrix,
[`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md)
builds the indicator matrix `B_t` for each treatment and the projector
`K_t` onto its null space. The vectors `b` that identify an effect are
the binary solutions of `b K_t = 0`.

``` r

KB <- projection_matrices(R, treatments)
b  <- solve_b_pairs(KB)

target_populations(b)
#> $cheap_mid
#> $cheap_mid$`A4+A7`
#> [1] "A4" "A7"
#> 
#> 
#> $cheap_dear
#> $cheap_dear$A5
#> [1] "A5"
#> 
#> 
#> $mid_dear
#> $mid_dear$`A6+A7`
#> [1] "A6" "A7"
```

All three pairwise contrasts are identified, but each for a
**different** population:

- `cheap` vs `mid` applies to teams whose adherence set contains both,
  that is `A4 = {cheap, mid}` or `A7 = {cheap, mid, dear}`.
- `cheap` vs `dear` applies only to teams whose set is *exactly*
  `A5 = {cheap, dear}`. Teams that would also consider `mid` are
  excluded, because for them a price change can divert choice to `mid`
  rather than between the two alternatives being compared.
- `mid` vs `dear` applies to teams containing both.

Reading this off is the point of the framework. An effect estimate means
little without knowing whose effect it is.

## Estimating

[`simulate_adherence()`](https://kazemia.github.io/catIV/reference/simulate_adherence.md)
generates data by exactly this mechanism, so we know the answer in
advance. Treatment response probabilities are 0.20, 0.35 and 0.50, and
the confounder shifts the outcome as well as the adherence set — so a
crude comparison will be biased.

``` r

set.seed(1)
d <- simulate_adherence(
  n = 40000, instrument = instrument, sets = sets, treatments = treatments,
  v_levels = c(0, 1), v_probs = c(0.5, 0.5),
  v_to_set = list("0" = c(3, 1, 1, 2, 4, 1, 3),
                  "1" = c(1, 3, 2, 1, 2, 3, 4)),
  v_on_y = c(0.10, 0.30), t_on_y = c(0.20, 0.35, 0.50)
)
head(d)
#>      Z     T Y  A
#> 1 2021   mid 1 A2
#> 2 2021   mid 1 A2
#> 3 2022 cheap 0 A1
#> 4 2020   mid 0 A7
#> 5 2019   mid 1 A6
#> 6 2020   mid 1 A4
```

Estimation is two observable quantities and then algebra. `P_Z` is the
probability of each treatment at each instrument value:

``` r

p_z <- estimate_p_z(d, "Z", "T", instrument_levels = rownames(R),
                    treatments = treatments)
round(p_z, 3)
#>      cheap   mid  dear
#> 2019 0.642 0.258 0.100
#> 2020 0.327 0.583 0.090
#> 2021 0.136 0.567 0.297
#> 2022 0.128 0.226 0.646
```

Treatment tracks the price ranking, which is the instrument doing its
work. `Q_Z` combines outcome and treatment, and
[`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md)
gives the size of each target population:

``` r

q_z     <- estimate_q_z(d, "Z", "T", "Y", instrument_levels = rownames(R),
                        treatments = treatments)
p_sigma <- estimate_p_sigma(p_z, KB, b)

lapply(p_sigma$average, round, 3)
#> $cheap_mid
#> A4+A7 
#> 0.316 
#> 
#> $cheap_dear
#>    A5 
#> 0.199 
#> 
#> $mid_dear
#> A6+A7 
#> 0.349
```

``` r

late <- estimate_late(q_z, KB, b, p_sigma)
lapply(late, round, 3)
#> $cheap_mid
#>  A4+A7 
#> -0.152 
#> 
#> $cheap_dear
#>     A5 
#> -0.308 
#> 
#> $mid_dear
#> A6+A7 
#> -0.15
```

True contrasts are `-0.15`, `-0.30` and `-0.15`. Compare with a crude
comparison of outcomes by treatment received:

``` r

crude <- tapply(d$Y, d$T, mean)
round(c(cheap_mid  = crude[["cheap"]] - crude[["mid"]],
        cheap_dear = crude[["cheap"]] - crude[["dear"]],
        mid_dear   = crude[["mid"]]   - crude[["dear"]]), 3)
#>  cheap_mid cheap_dear   mid_dear 
#>     -0.210     -0.342     -0.132
```

The crude estimates are pulled toward zero because teams inclined to the
dearer options also have better outcomes for reasons unrelated to
treatment.

## Why not just compare two at a time?

A natural shortcut is to throw away everyone who took a third
alternative, collapse the instrument to “is A cheaper than B?”, and run
the usual Wald ratio.
[`naive_iv()`](https://kazemia.github.io/catIV/reference/naive_iv.md)
does exactly that:

``` r

round(c(
  cheap_mid  = naive_iv("cheap_mid",  d, instrument, "Z", "T", "Y"),
  cheap_dear = naive_iv("cheap_dear", d, instrument, "Z", "T", "Y"),
  mid_dear   = naive_iv("mid_dear",   d, instrument, "Z", "T", "Y")
), 3)
#>  cheap_mid cheap_dear   mid_dear 
#>     -0.159     -0.320     -0.126
```

These are much better than the crude comparison, because the instrument
still removes the confounding. What they do not fix is *which* movement
is being measured. Making `cheap` cheaper does not only move teams
between `cheap` and `mid`; it also pulls in teams who would otherwise
have taken `dear`. The two-at-a-time analysis attributes that movement
to the `cheap` versus `mid` contrast, so the estimate is a blend of
contrasts rather than the one asked for.

How much that matters depends on the design. Here the effects are
homogeneous across decision teams, so blending contrasts is mild and the
naive estimates land within a few hundredths — the `mid` versus `dear`
estimate is the most visibly off. With effects that differ across
alternatives or across teams the discrepancy grows, and in either case
the naive estimate has no clean interpretation: there is no population
for which it is the average effect.

## Restricting the adherence sets

Some adherence sets may be clinically implausible, and ruling them out
is an extra assumption. It changes what is identified and for whom — not
always helpfully. Here is the effect of dropping each set in turn:

``` r

describe <- function(keep) {
  tp <- target_populations(
    solve_b_pairs(projection_matrices(response_matrix(sets[keep], instrument),
                                      treatments))
  )
  vapply(tp, function(x) if (length(x) == 0) "none" else
    paste(names(x), collapse = " or "), character(1))
}

t(vapply(names(sets), function(drop) describe(setdiff(names(sets), drop)),
         character(3)))
#>    cheap_mid cheap_dear mid_dear
#> A1 "A4+A7"   "A5"       "A6+A7" 
#> A2 "A4+A7"   "A5"       "A6+A7" 
#> A3 "A4+A7"   "A5"       "A6+A7" 
#> A4 "A7"      "A5"       "A6+A7" 
#> A5 "A4+A7"   "none"     "A6+A7" 
#> A6 "A4+A7"   "A5"       "A7"    
#> A7 "A4"      "A5"       "A6"
```

Dropping `A5` removes the `cheap` vs `dear` contrast altogether: `A5`
was the only population for which it was identified. Dropping `A4`
narrows `cheap` vs `mid` from `A4 or A7` to `A7` alone. Restriction here
costs information rather than adding it, because with three treatments
and four instrument values all three contrasts are identified already.

The benefit appears with more alternatives, and it is **descriptive**.
With five treatments there are 31 adherence sets, and the target
population of an identified contrast can be a union of many of them —
arithmetically fine, impossible to describe to a clinician:

``` r

# The design in paper 2: five TNF inhibitors, ten procurement periods.
target_populations(solve_b_pairs(projection_matrices(R_full, TNFi)))$Cer_Inf
#> $`A12+A18+A23+A24+A27+A28+A30+A31`
#> [1] "A12" "A18" "A23" "A24" "A27" "A28" "A30" "A31"
```

Restricting to the 14 adherence sets clinicians considered realistic
leaves the same seven contrasts identified, but the target population
for that contrast becomes a single set — decision teams open to any of
the five drugs. A describable population is what makes the estimate
usable, which is why the restriction was worth the assumption.

## Next steps

- [`vignette("target-trials")`](https://kazemia.github.io/catIV/articles/target-trials.md)
  covers covariate adjustment, bootstrap confidence intervals, and
  describing the target population as a trial.
- [`vignette("heterogeneous-effects")`](https://kazemia.github.io/catIV/articles/heterogeneous-effects.md)
  covers effects that vary with a covariate.
