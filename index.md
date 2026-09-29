# catIV

Instrumental variable analysis when the instrument is **ordinal** and
the treatment is **categorical** — more than two alternatives, ranked by
how much the instrument encourages each one.

Standard IV methods handle a binary treatment and identify an effect for
compliers. With three or more alternatives, “complier” is no longer a
single group: a decision team might be willing to consider alternatives
A and B but never C, and the instrument moves them between A and B only.
`catIV` formalises this through **adherence sets** — the set of
alternatives a decision team would actually consider — works out which
contrasts are identifiable given the instrument, and says exactly which
subpopulation each identified effect applies to.

The methods come from three papers on comparing TNF-α inhibitors for
rheumatoid arthritis, using national drug-price rankings as the
instrument.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("kazemia/catIV")
```

## A worked example

Three treatments ranked by price, with the ranking changing across four
procurement periods. The ranking is the instrument.

``` r

library(catIV)

treatments <- c("cheap", "mid", "dear")
instrument <- list(
  "2019" = c("cheap", "mid",   "dear"),
  "2020" = c("mid",   "cheap", "dear"),
  "2021" = c("mid",   "dear",  "cheap"),
  "2022" = c("dear",  "mid",   "cheap")
)
```

### 1. What is identifiable?

With three treatments there are seven possible adherence sets. Which
contrasts can this instrument identify, and for whom?

``` r

sets <- adherence_sets(treatments)     # A1 = {cheap} ... A7 = {cheap, mid, dear}
R    <- response_matrix(sets, instrument)
KB   <- projection_matrices(R, treatments)
b    <- solve_b_pairs(KB)

target_populations(b)
#> $cheap_mid
#> $cheap_mid$`A4+A7`
#> [1] "A4" "A7"
#>
#> $cheap_dear
#> $cheap_dear$A5
#> [1] "A5"
#>
#> $mid_dear
#> $mid_dear$`A6+A7`
#> [1] "A6" "A7"
```

All three pairwise contrasts are identified, each for a different
population. `cheap` versus `mid` applies to decision teams considering
**both** of them (`A4 = {cheap, mid}` or `A7 = {cheap, mid, dear}`),
whereas `cheap` versus `dear` applies only to those considering
**exactly** `{cheap, dear}` and not `mid`. That is the categorical
analogue of “compliers”, and reading it off is the point of the
framework.

### 2. Estimate

``` r

set.seed(1)
d <- simulate_adherence(
  n = 5000, instrument = instrument, sets = sets, treatments = treatments,
  v_levels = c(0, 1), v_probs = c(0.5, 0.5),
  v_to_set = list("0" = c(3, 1, 1, 2, 4, 1, 3), "1" = c(1, 3, 2, 1, 2, 3, 4)),
  v_on_y = c(0.10, 0.30), t_on_y = c(0.20, 0.35, 0.50)
)

p_z <- estimate_p_z(d, "Z", "T", instrument_levels = rownames(R),
                    treatments = treatments)
round(p_z, 3)
#>      cheap   mid  dear
#> 2019 0.631 0.259 0.110
#> 2020 0.328 0.582 0.090
#> 2021 0.131 0.595 0.274
#> 2022 0.129 0.224 0.647
```

Treatment follows the price ranking, which is what makes the instrument
work.

``` r

q_z     <- estimate_q_z(d, "Z", "T", "Y", instrument_levels = rownames(R),
                        treatments = treatments)
p_sigma <- estimate_p_sigma(p_z, KB, b)

estimate_late(q_z, KB, b, p_sigma)
#> $cheap_mid
#>  A4+A7
#> -0.159
#>
#> $cheap_dear
#>     A5
#> -0.249
#>
#> $mid_dear
#>  A6+A7
#> -0.123
```

The simulated response probabilities were 0.20, 0.35 and 0.50, so the
true contrasts are −0.15, −0.30 and −0.15.

[`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md)
also tells you how large each target population is, which is how you
judge whether an estimate is worth reporting:

``` r

lapply(p_sigma$average, round, 3)
#> $cheap_mid
#> A4+A7
#> 0.318
#>
#> $cheap_dear
#>    A5
#> 0.186
#>
#> $mid_dear
#> A6+A7
#> 0.368
```

### 3. Confidence intervals

``` r

bootstrap_late(200, d, "Z", "T", "Y", KB, b, n_cores = 2)$ci
#>                  lower  upper n_kept
#> cheap_mid.A4+A7 -0.274 -0.055    200
#> cheap_dear.A5   -0.366 -0.121    200
#> mid_dear.A6+A7  -0.217 -0.031    200
```

## Other things the package does

- **Covariate adjustment.**
  `estimate_p_z(..., parametric = TRUE, covariates = ...)` standardises
  over observed covariates via multinomial regression, and
  [`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md)
  does the same for the outcome.
- **Target trial descriptions.**
  [`pseudo_population()`](https://kazemia.github.io/catIV/reference/pseudo_population.md)
  returns a weighted dataset on which an ordinary weighted comparison
  reproduces the IV estimate, so weighted baseline tables describe who
  the effect applies to.
- **Heterogeneous effects.**
  [`hiv_cate()`](https://kazemia.github.io/catIV/reference/hiv_cate.md)
  combines instrumental and covariate information to estimate
  conditional average treatment effects, keeping most of the efficiency
  of a covariate-adjusted model without its confounding bias.
- **Simulation.**
  [`simulate_adherence()`](https://kazemia.github.io/catIV/reference/simulate_adherence.md)
  generates data satisfying the identification assumptions exactly;
  [`simulate_categorical()`](https://kazemia.github.io/catIV/reference/simulate_categorical.md)
  generates data that only approximately satisfies them, for testing
  robustness.

## Function names used in the papers

The published code for the three papers used different names.

| Paper code | catIV | Notes |
|----|----|----|
| `GenerateA()` | [`adherence_sets()`](https://kazemia.github.io/catIV/reference/adherence_sets.md) |  |
| `T_decider()` | [`choose_treatment()`](https://kazemia.github.io/catIV/reference/choose_treatment.md) |  |
| `MakeR()` | [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md) |  |
| `MakeKB()` | [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md) | list elements `B_t`, `B_t_i` are now `B`, `B_plus` |
| `KbSolver()` | [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md) | papers 1 and 2, indexed by treatment pairs |
| `KbSolver()` | [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md) | paper 3, indexed by single treatments |
| `PiIdentifier()` | [`target_populations()`](https://kazemia.github.io/catIV/reference/target_populations.md) |  |
| `MakeP_Z()` | [`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md) |  |
| `MakeQ_Z()` | [`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md) |  |
| `MakeV_Z()` | [`estimate_v_z()`](https://kazemia.github.io/catIV/reference/estimate_v_z.md) |  |
| `P_SigmaIdentifier()` | [`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md) | `E1`, `E2`, `WA` are now `arm1`, `arm2`, `average` |
| `P_SigmaIdentifier()` | [`estimate_p_sigma_treatments()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma_treatments.md) | the paper 3 variant |
| `LATEIdentifier()` | [`estimate_late()`](https://kazemia.github.io/catIV/reference/estimate_late.md) | `RR`, `AverageProb` are now `scale`, `denominator` |
| `LATOIdentifier()` | [`estimate_latr()`](https://kazemia.github.io/catIV/reference/estimate_latr.md) |  |
| `BSCICalculator()` | [`bootstrap_late()`](https://kazemia.github.io/catIV/reference/bootstrap_late.md) | `Data.complete = FALSE` is now just passing a list |
| `PseudoPopulator()` | [`pseudo_population()`](https://kazemia.github.io/catIV/reference/pseudo_population.md) | `Pi_index` is now `solution`, given by name |
| `NaiveIV()` | [`naive_iv()`](https://kazemia.github.io/catIV/reference/naive_iv.md) |  |
| `Pi_prob_estimator()` | [`conditional_p_sigma()`](https://kazemia.github.io/catIV/reference/conditional_p_sigma.md) |  |
| `CIV_estimator()` | [`conditional_latr()`](https://kazemia.github.io/catIV/reference/conditional_latr.md) |  |
| `Confounded_estimator()` | [`confounded_response()`](https://kazemia.github.io/catIV/reference/confounded_response.md) |  |
| `BinarySimulator()` | [`simulate_binary()`](https://kazemia.github.io/catIV/reference/simulate_binary.md) | `OR` is now `logit` |
| `CatSimulator()` (paper 1) | [`simulate_adherence()`](https://kazemia.github.io/catIV/reference/simulate_adherence.md) |  |
| `CatSimulator2()` (paper 1) | [`simulate_categorical()`](https://kazemia.github.io/catIV/reference/simulate_categorical.md) | paper 2’s `CatSimulator()` is the same with `intercept = 0` |
| `C_probability()` and friends | [`binary_iv_truth()`](https://kazemia.github.io/catIV/reference/binary_iv_truth.md) | all four quantities from one marginalisation |

## Differences from the published code

Results are numerically identical to the paper code, verified against it
across the designs in all three papers. A few conventions changed:

- `P_Z` and `Q_Z` are matrices with instrument values as row names,
  rather than data frames whose first column holds the instrument.
  Nothing depends on column position.
- The solvers always return a matrix; a contrast with no identifiable
  effect gives a matrix with zero rows.
- Solutions are named by the adherence sets they select, such as
  `"A4+A7"`. Solution *order* therefore differs from the paper code, so
  index-based references like `Pi_index = 1` do not carry over.
- Every instrument value must be observed in the data. The paper code
  silently produced a shorter `P_Z`.
- [`hiv_cate()`](https://kazemia.github.io/catIV/reference/hiv_cate.md)
  weights the bias regression by the covariate density to the first
  power, as specified in paper 3. The analysis script used the density
  squared; pass `density_power = 2` to reproduce those published
  numbers.
- [`binary_iv_truth()`](https://kazemia.github.io/catIV/reference/binary_iv_truth.md)
  marginalises every quantity over both confounders. The paper 3 helper
  `T_probability()` omitted the observed confounder, which `v_on_t = 0`
  reproduces.

## Licence

MIT
