hte_fit <- function(n = 2500, seed = 42) {
  treatments <- c("t1", "t2", "t3")
  instrument <- list(`1` = c("t1", "t2", "t3"), `2` = c("t2", "t1", "t3"),
                     `3` = c("t2", "t3", "t1"), `4` = c("t3", "t2", "t1"))
  R <- response_matrix(adherence_sets(treatments), instrument)
  KB <- projection_matrices(R, treatments)
  b <- solve_b_treatments(KB)
  sets <- adherence_sets(treatments)

  set.seed(seed)
  X <- stats::runif(n, 20, 80)
  U <- stats::rbinom(n, 1, 0.5)
  w <- outer(seq_along(sets), seq_len(n),
             function(i, j) exp(0.25 * i + 0.03 * X[j] * (i %% 3) +
                                  0.8 * U[j] * (i %% 2)))
  a <- apply(w, 2, function(p) sample(seq_along(sets), 1, prob = p))
  z <- sample(seq_along(instrument), n, TRUE)
  trt <- vapply(seq_len(n),
                function(i) choose_treatment(instrument[[z[i]]], sets[[a[i]]]),
                character(1))
  base <- c(t1 = 0.10, t2 = 0.35, t3 = 0.60)
  Y <- 0.2 + base[trt] + 0.005 * X + 0.9 * U + stats::rnorm(n, 0, 0.3)
  d <- data.frame(Z = factor(z, levels = seq_along(instrument)),
                  T = factor(trt, levels = treatments),
                  X = X, Y = as.numeric(Y))

  list(treatments = treatments, instrument = instrument,
       z_levels = as.character(seq_along(instrument)),
       KB = KB, b = b, d = d,
       den = nnet::multinom(T ~ Z + X, data = d, trace = FALSE),
       num = stats::lm(Y ~ Z * T * X, data = d),
       conf = stats::lm(Y ~ T * X, data = d))
}

test_that("conditional_p_z gives a row-stochastic matrix of the right shape", {
  f <- hte_fit()
  p <- conditional_p_z(50, f$den, f$z_levels, f$treatments)
  expect_equal(dim(p), c(length(f$z_levels), length(f$treatments)))
  expect_identical(rownames(p), f$z_levels)
  expect_identical(colnames(p), f$treatments)
  expect_equal(unname(rowSums(p)), rep(1, length(f$z_levels)))
  expect_true(all(p >= 0))
})

test_that("conditional_p_z varies with the covariate", {
  f <- hte_fit()
  expect_false(isTRUE(all.equal(conditional_p_z(25, f$den, f$z_levels, f$treatments),
                                conditional_p_z(75, f$den, f$z_levels, f$treatments))))
})

test_that("conditional_p_sigma matches the unconditional form of the algebra", {
  f <- hte_fit()
  p_z <- conditional_p_z(50, f$den, f$z_levels, f$treatments)
  expect_equal(conditional_p_sigma(50, f$den, f$KB, f$b, f$z_levels, f$treatments),
               estimate_p_sigma_treatments(p_z, f$KB, f$b))
})

test_that("conditional_p_sigma stays a probability across the covariate range", {
  f <- hte_fit()
  for (x in c(25, 40, 55, 70)) {
    ps <- conditional_p_sigma(x, f$den, f$KB, f$b, f$z_levels, f$treatments)
    vals <- unlist(ps)
    expect_true(all(vals > -1e-6 & vals < 1 + 1e-6))
  }
})

test_that("conditional_latr returns one response per treatment and solution", {
  f <- hte_fit()
  out <- conditional_latr(50, f$num, f$den, f$KB, f$b, f$z_levels, f$treatments)
  expect_named(out, f$treatments)
  for (t in f$treatments) {
    expect_identical(names(out[[t]]), rownames(f$b[[t]]))
  }
})

test_that("hiv_cate returns a debiased response per covariate value and treatment", {
  f <- hte_fit()
  xs <- seq(30, 70, length.out = 6)
  dens <- stats::density(f$d$X, from = min(xs), to = max(xs), n = length(xs))$y
  fit <- hiv_cate(xs, f$num, f$den, f$conf, f$KB, f$b, density = dens,
                  instrument_levels = f$z_levels, treatments = f$treatments)
  expect_equal(nrow(fit), length(xs) * length(f$treatments))
  expect_named(fit, c("x", "treatment", "confounded", "bias", "response"))
  expect_equal(fit$response, fit$confounded - fit$bias)
  expect_s3_class(attr(fit, "bias_model"), "lm")
  expect_true(all(c("p_sigma", "latr", "w") %in% names(attr(fit, "bias_data"))))
})

test_that("hiv_cate weights use the density raised to density_power", {
  f <- hte_fit()
  xs <- seq(30, 70, length.out = 5)
  dens <- stats::density(f$d$X, from = min(xs), to = max(xs), n = length(xs))$y
  for (power in c(1, 2)) {
    bd <- attr(hiv_cate(xs, f$num, f$den, f$conf, f$KB, f$b, density = dens,
                        instrument_levels = f$z_levels,
                        treatments = f$treatments, density_power = power),
               "bias_data")
    expect_equal(bd$w, bd$density^power * bd$p_sigma^2)
  }
})

test_that("density must line up with the covariate values", {
  f <- hte_fit()
  expect_error(
    hiv_cate(c(30, 50), f$num, f$den, f$conf, f$KB, f$b, density = 1,
             instrument_levels = f$z_levels, treatments = f$treatments),
    "one value per element"
  )
})

test_that("hiv_contrast differences two treatments and average_contrast collapses it", {
  f <- hte_fit()
  xs <- seq(30, 70, length.out = 6)
  dens <- stats::density(f$d$X, from = min(xs), to = max(xs), n = length(xs))$y
  fit <- hiv_cate(xs, f$num, f$den, f$conf, f$KB, f$b, density = dens,
                  instrument_levels = f$z_levels, treatments = f$treatments)
  con <- hiv_contrast(fit, "t1", "t3")
  expect_equal(con$estimate,
               fit$response[fit$treatment == "t1"] -
                 fit$response[fit$treatment == "t3"])
  expect_equal(hiv_contrast(fit, "t3", "t1")$estimate, -con$estimate)
  expect_equal(average_contrast(con, dens),
               sum(con$estimate * dens) / sum(dens))
  expect_error(hiv_contrast(fit, "t1", "nope"), "Both treatments")
})
