# Effects that vary with a covariate

``` r

library(catIV)
```

## Two estimators, each flawed

Suppose the treatment effect varies with age, and age also influences
which treatments a decision team will consider. Two obvious approaches
both fail:

- **Adjust for observed covariates.** Efficient, and it gives an
  estimate at every age, but it is biased by whatever confounding you
  did not observe.
- **Use the instrument, conditioning on age.** Unbiased, but the
  denominator is the probability of being in the target population *at
  that age*. Where few decision teams are open to both alternatives,
  that probability approaches zero and the estimate becomes unusable.

Paper 3’s estimator combines them. It models the *difference* between
the two — the bias of the covariate-adjusted estimator — and subtracts
the fitted bias from the efficient estimate. Where the instrument is
informative the bias is pinned down well; where it is not, the weighting
lets neighbouring regions carry the estimate.

## Setup

``` r

treatments <- c("t1", "t2", "t3")
instrument <- list(`1` = c("t1", "t2", "t3"), `2` = c("t2", "t1", "t3"),
                   `3` = c("t2", "t3", "t1"), `4` = c("t3", "t2", "t1"))

sets <- adherence_sets(treatments)
R    <- response_matrix(sets, instrument)
KB   <- projection_matrices(R, treatments)
```

For conditional quantities the relevant solutions are the
**single-treatment** ones, solving `b K_t = 0` for one treatment at a
time. These identify local average treatment *responses*,
`E(Y^{T=t} | b[A] = 1, X)`, which the estimator works with one arm at a
time before differencing.

``` r

b <- solve_b_treatments(KB)
lengths(lapply(b, rownames))
#> t1 t2 t3 
#>  7  5  7
```

Data where the adherence set depends on age and on an unobserved
confounder, and the effect of each treatment varies linearly with age:

``` r

set.seed(314)
n <- 4000
X <- runif(n, 20, 80)
U <- rbinom(n, 1, 0.5)          # unobserved
V <- rbinom(n, 1, 0.5)          # observed

w <- outer(seq_along(sets), seq_len(n),
           function(i, j) exp(0.25 * i + 0.03 * X[j] * (i %% 3) +
                                0.8 * U[j] * (i %% 2) +
                                0.5 * V[j] * (i %% 3 == 0)))
a <- apply(w, 2, function(p) sample(seq_along(sets), 1, prob = p))
z <- sample(seq_along(instrument), n, TRUE)
trt <- vapply(seq_len(n),
              function(i) choose_treatment(instrument[[z[i]]], sets[[a[i]]]),
              character(1))

base  <- c(t1 = 0.10, t2 = 0.35, t3 = 0.60)
slope <- c(t1 = 0.30, t2 = 0.10, t3 = -0.20)
Y <- 0.2 + base[trt] + slope[trt] * (X - 50) / 50 +
  0.005 * X + 0.9 * U + 0.4 * V + rnorm(n, 0, 0.3)

d <- data.frame(Z = factor(z, levels = seq_along(instrument)),
                T = factor(trt, levels = treatments),
                X = X, Y = as.numeric(Y), V = V)
```

## Where is the instrument informative?

Before estimating anything, check whether the target population exists
across the range of ages.
[`conditional_p_sigma()`](https://kazemia.github.io/catIV/reference/conditional_p_sigma.md)
gives its probability at any value of the covariate:

``` r

den <- nnet::multinom(T ~ Z + X, data = d, trace = FALSE)

ages <- seq(25, 75, length.out = 11)
sizes <- t(vapply(ages, function(x) {
  ps <- conditional_p_sigma(x, den, KB, b, levels(d$Z), treatments)
  vapply(treatments, function(t) max(ps[[t]]), numeric(1))
}, numeric(length(treatments))))
rownames(sizes) <- ages
round(sizes, 3)
#>       t1    t2    t3
#> 25 0.723 0.503 0.700
#> 30 0.734 0.490 0.707
#> 35 0.745 0.477 0.714
#> 40 0.756 0.464 0.720
#> 45 0.766 0.451 0.726
#> 50 0.776 0.439 0.733
#> 55 0.785 0.426 0.739
#> 60 0.795 0.413 0.744
#> 65 0.803 0.401 0.750
#> 70 0.812 0.388 0.755
#> 75 0.820 0.376 0.760
```

Where a column approaches zero, the IV estimator for that treatment has
a vanishing denominator. That is the positivity problem the combined
estimator exists to work around.

## The estimator

Three models: the treatment given instrument and covariate, the outcome
given all three, and a covariate-adjusted outcome model that does *not*
use the instrument.

``` r

num  <- lm(Y ~ Z * T * X, data = d)
conf <- lm(Y ~ T * X + V, data = d)
```

The estimator also needs the density of the covariate, because the
weights give more influence to regions where the data are dense.

``` r

dens <- density(d$X, from = min(ages), to = max(ages), n = length(ages))$y

fit <- hiv_cate(ages, num, den, conf, KB, b, density = dens,
                instrument_levels = levels(d$Z), treatments = treatments)
head(fit)
#>    x treatment confounded       bias  response
#> 1 25        t1  0.7375143 -0.2074113 0.9449256
#> 2 30        t1  0.7955205 -0.2030171 0.9985376
#> 3 35        t1  0.8535267 -0.1988865 1.0524132
#> 4 40        t1  0.9115329 -0.1950111 1.1065439
#> 5 45        t1  0.9695390 -0.1913827 1.1609217
#> 6 50        t1  1.0275452 -0.1879931 1.2155384
```

`response` is the debiased conditional average treatment response: the
covariate-adjusted estimate minus the fitted bias. Contrasts are
differences between two of them:

``` r

con <- hiv_contrast(fit, "t1", "t3")
truth <- (base[["t1"]] - base[["t3"]]) +
  (slope[["t1"]] - slope[["t3"]]) * (ages - 50) / 50

round(data.frame(age = ages, estimate = con$estimate, truth = truth), 3)
#>    age estimate truth
#> 1   25   -0.767 -0.75
#> 2   30   -0.714 -0.70
#> 3   35   -0.661 -0.65
#> 4   40   -0.607 -0.60
#> 5   45   -0.554 -0.55
#> 6   50   -0.501 -0.50
#> 7   55   -0.448 -0.45
#> 8   60   -0.395 -0.40
#> 9   65   -0.342 -0.35
#> 10  70   -0.289 -0.30
#> 11  75   -0.235 -0.25
```

And
[`average_contrast()`](https://kazemia.github.io/catIV/reference/average_contrast.md)
collapses it to a single number, averaging over the covariate
distribution:

``` r

round(average_contrast(con, dens), 3)
#> [1] -0.502
round(sum(truth * dens) / sum(dens), 3)
#> [1] -0.501
```

## The weighting

The bias regression weights each observation by
`density^density_power * p_sigma^2`, so regions where the instrument is
informative and the data are dense carry more of the fit.

**`density_power` defaults to 1**, matching the variance argument in
paper 3, which makes the variance proportional to
`1 / (p_sigma^2 * density)`. The analysis code used to produce the
published clinical results squared the density instead. To reproduce
those numbers exactly, set `density_power = 2`:

``` r

published <- hiv_cate(ages, num, den, conf, KB, b, density = dens,
                      instrument_levels = levels(d$Z),
                      treatments = treatments, density_power = 2)

round(hiv_contrast(published, "t1", "t3")$estimate - con$estimate, 4)
#>  [1] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e-04 1e-04 1e-04 1e-04 1e-04
```

The fitted bias model and the assembled per-solution data are attached
as attributes, which is where to look when a fit is behaving oddly:

``` r

head(attr(fit, "bias_data")[, c("x", "treatment", "solution",
                                "p_sigma", "latr", "confounded", "w")], 4)
#>    x treatment solution   p_sigma      latr confounded           w
#> 1 25        t1       A5 0.4221202 0.9240165  0.7375143 0.002869295
#> 2 30        t1       A5 0.4354333 0.9832159  0.7955205 0.003066018
#> 3 35        t1       A5 0.4487462 1.0424275  0.8535267 0.003365037
#> 4 40        t1       A5 0.4620332 1.1016524  0.9115329 0.003507163
```

## Inference

There is no closed form. Inference means bootstrapping the whole
procedure — refitting all three models and re-running the estimator on
each resample — and taking pointwise quantiles across the covariate
range. That is expensive, which is the practical cost of the method.
