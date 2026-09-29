# Package index

## Describing the design

Enumerate the adherence sets a decision team could hold, and record
which treatment each one receives under each value of the instrument.

- [`adherence_sets()`](https://kazemia.github.io/catIV/reference/adherence_sets.md)
  : Enumerate all possible adherence sets
- [`choose_treatment()`](https://kazemia.github.io/catIV/reference/choose_treatment.md)
  : The choice function implied by categorical monotonicity
- [`response_matrix()`](https://kazemia.github.io/catIV/reference/response_matrix.md)
  : Build the response matrix

## Identification

Work out which contrasts the instrument can identify, and which
subpopulation each identified effect applies to.

- [`projection_matrices()`](https://kazemia.github.io/catIV/reference/projection_matrices.md)
  : Indicator matrices and null-space projectors
- [`solve_b_pairs()`](https://kazemia.github.io/catIV/reference/solve_b_pairs.md)
  : Identify the target populations of local average treatment effects
- [`solve_b_treatments()`](https://kazemia.github.io/catIV/reference/solve_b_treatments.md)
  : Identify the target populations of local average treatment responses
- [`target_populations()`](https://kazemia.github.io/catIV/reference/target_populations.md)
  : Describe the identified target populations

## Estimation

Estimate the observable quantities from data, with or without covariate
adjustment, and combine them into local effects.

- [`estimate_p_z()`](https://kazemia.github.io/catIV/reference/estimate_p_z.md)
  : Estimate the conditional treatment probabilities P_Z
- [`estimate_q_z()`](https://kazemia.github.io/catIV/reference/estimate_q_z.md)
  : Estimate the outcome-probability products Q_Z
- [`estimate_v_z()`](https://kazemia.github.io/catIV/reference/estimate_v_z.md)
  : Estimate the conditional outcome variances V_Z
- [`estimate_p_sigma()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma.md)
  : Probability of belonging to each identified target population
- [`estimate_p_sigma_treatments()`](https://kazemia.github.io/catIV/reference/estimate_p_sigma_treatments.md)
  : Probability of belonging to each single-treatment target population
- [`estimate_late()`](https://kazemia.github.io/catIV/reference/estimate_late.md)
  : Identify local average treatment effects
- [`estimate_latr()`](https://kazemia.github.io/catIV/reference/estimate_latr.md)
  : Identify local average treatment responses
- [`naive_iv()`](https://kazemia.github.io/catIV/reference/naive_iv.md)
  : A naive instrumental variable estimate ignoring the other treatments

## Inference and interpretation

Bootstrap confidence intervals, and the weighted pseudo-population that
describes the target trial an estimate emulates.

- [`bootstrap_late()`](https://kazemia.github.io/catIV/reference/bootstrap_late.md)
  : Bootstrap confidence intervals for local average treatment effects
- [`pseudo_population()`](https://kazemia.github.io/catIV/reference/pseudo_population.md)
  : Build the pseudo-population behind an identified effect

## Heterogeneous treatment effects

Conditional versions of the estimators, and the estimator combining
instrumental with covariate information.

- [`conditional_p_z()`](https://kazemia.github.io/catIV/reference/conditional_p_z.md)
  : Conditional treatment probabilities from a fitted model
- [`conditional_p_sigma()`](https://kazemia.github.io/catIV/reference/conditional_p_sigma.md)
  : Conditional probability of belonging to each target population
- [`conditional_latr()`](https://kazemia.github.io/catIV/reference/conditional_latr.md)
  : Conditional local average treatment responses
- [`confounded_response()`](https://kazemia.github.io/catIV/reference/confounded_response.md)
  : Conditional treatment responses from a confounded outcome model
- [`hiv_cate()`](https://kazemia.github.io/catIV/reference/hiv_cate.md)
  : Heterogeneous instrumental variable estimator of conditional
  responses
- [`hiv_contrast()`](https://kazemia.github.io/catIV/reference/hiv_contrast.md)
  : Contrast two treatments from a heterogeneous IV fit
- [`average_contrast()`](https://kazemia.github.io/catIV/reference/average_contrast.md)
  : Average a conditional effect over the covariate distribution

## Simulation

Data-generating processes for checking estimators against a known truth.

- [`simulate_adherence()`](https://kazemia.github.io/catIV/reference/simulate_adherence.md)
  : Simulate from the adherence set model
- [`simulate_categorical()`](https://kazemia.github.io/catIV/reference/simulate_categorical.md)
  : Simulate a categorical treatment chosen by a multinomial model
- [`simulate_binary()`](https://kazemia.github.io/catIV/reference/simulate_binary.md)
  : Simulate a binary treatment with a binary instrument
- [`binary_iv_truth()`](https://kazemia.github.io/catIV/reference/binary_iv_truth.md)
  : True principal stratum probabilities for a binary instrument
