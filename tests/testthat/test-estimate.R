toy <- function() {
  # Two instrument values, two treatments, hand-countable.
  data.frame(
    z   = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
    trt = c("a", "a", "a", "b", "a", "b", "b", "b", "b", "b"),
    y   = c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1)
  )
}

test_that("estimate_p_z returns row-stochastic cell proportions", {
  p <- estimate_p_z(toy(), "z", "trt")
  expect_equal(dim(p), c(2L, 2L))
  expect_identical(rownames(p), c("1", "2"))
  expect_identical(colnames(p), c("a", "b"))
  expect_equal(unname(rowSums(p)), c(1, 1))
  expect_equal(unname(p["1", "a"]), 3 / 4)
  expect_equal(unname(p["2", "b"]), 5 / 6)
})

test_that("estimate_q_z equals within-cell outcome sums over the instrument count", {
  q <- estimate_q_z(toy(), "z", "trt", "y")
  # Z = 1: outcomes for a are 1, 0, 1 (sum 2); n at Z = 1 is 4.
  expect_equal(unname(q["1", "a"]), 2 / 4)
  expect_equal(unname(q["1", "b"]), 1 / 4)
  # Z = 2: outcomes for b are 1, 1, 0, 1, 1 (sum 4); n at Z = 2 is 6.
  expect_equal(unname(q["2", "b"]), 4 / 6)
})

test_that("Q_Z divided by P_Z recovers the conditional outcome mean", {
  d <- toy()
  p <- estimate_p_z(d, "z", "trt")
  q <- estimate_q_z(d, "z", "trt", "y")
  expect_equal(unname(q["1", "a"] / p["1", "a"]),
               mean(d$y[d$z == 1 & d$trt == "a"]))
  expect_equal(unname(q["2", "b"] / p["2", "b"]),
               mean(d$y[d$z == 2 & d$trt == "b"]))
})

test_that("empty cells give zero rather than NA", {
  d <- toy()
  d <- d[!(d$z == 1 & d$trt == "b"), ]
  expect_equal(unname(estimate_p_z(d, "z", "trt")["1", "b"]), 0)
  expect_equal(unname(estimate_q_z(d, "z", "trt", "y")["1", "b"]), 0)
  expect_equal(unname(estimate_v_z(d, "z", "trt", "y")["1", "b"]), 0)
})

test_that("estimate_v_z scales the within-cell variance by the squared share", {
  d <- toy()
  v <- estimate_v_z(d, "z", "trt", "y")
  y11 <- d$y[d$z == 1 & d$trt == "a"]
  expect_equal(unname(v["1", "a"]), var(y11) * (3 / 4)^2)
  # A cell with one observation has no estimable variance and contributes zero.
  expect_equal(unname(v["1", "b"]), 0)
})

test_that("covariate adjustment is refused rather than silently ignored", {
  d <- toy(); d$c <- rnorm(nrow(d))
  expect_error(estimate_p_z(d, "z", "trt", covariates = "c"), "parametric")
  expect_error(estimate_q_z(d, "z", "trt", "y", covariates = "c"), "parametric")
  expect_error(estimate_q_z(d, "z", "trt", "y", parametric = TRUE), "`p_z`")
})

test_that("an unobserved instrument value is an error, not a short matrix", {
  expect_error(
    estimate_p_z(toy(), "z", "trt", instrument_levels = c("1", "2", "3")),
    "No observations for instrument value"
  )
  expect_error(estimate_p_z(toy(), "nope", "trt"), "not present")
})

# --- identification-side estimators ----------------------------------------

# Generate data the way the model assumes: draw an adherence set, draw an
# instrument value, then let categorical monotonicity pick the treatment. The
# outcome depends on the treatment and, separately, on the adherence set, so
# the adherence set confounds the crude comparison while every treatment
# effect stays homogeneous. Each identified LATE must therefore equal the true
# contrast, whatever population it is local to.
thesis_fit <- function(n = 30000, seed = 1) {
  treatments <- c("t1", "t2", "t3")
  instrument <- list(`1` = c("t1", "t2", "t3"), `2` = c("t2", "t1", "t3"),
                     `3` = c("t2", "t3", "t1"), `4` = c("t3", "t2", "t1"))
  R <- response_matrix(adherence_sets(treatments), instrument)
  KB <- projection_matrices(R, treatments)
  sets <- adherence_sets(treatments)
  truth <- c(t1 = 0.30, t2 = 0.50, t3 = 0.70)
  # Adherence sets that never respond to the instrument get little weight;
  # A5 and A7 carry the identified populations.
  set_weight <- c(A1 = 1, A2 = 1, A3 = 1, A4 = 2, A5 = 3, A6 = 2, A7 = 4)
  set_shift  <- c(A1 = -0.10, A2 = 0.05, A3 = 0.10, A4 = -0.05,
                  A5 = 0.08, A6 = -0.08, A7 = 0.00)

  set.seed(seed)
  a <- sample(names(sets), n, TRUE, prob = set_weight)
  z <- sample(seq_along(instrument), n, TRUE)
  trt <- vapply(seq_len(n),
                function(i) choose_treatment(instrument[[z[i]]], sets[[a[i]]]),
                character(1))
  d <- data.frame(z = z, trt = trt, a = a,
                  y = stats::rbinom(n, 1, truth[trt] + set_shift[a]))
  list(treatments = treatments, R = R, KB = KB, d = d, truth = truth,
       b = solve_b_pairs(KB),
       p_z = estimate_p_z(d, "z", "trt", instrument_levels = rownames(R),
                          treatments = treatments),
       q_z = estimate_q_z(d, "z", "trt", "y", instrument_levels = rownames(R),
                          treatments = treatments))
}

test_that("estimate_p_sigma averages the two arms by their marginal shares", {
  f <- thesis_fit()
  ps <- estimate_p_sigma(f$p_z, f$KB, f$b)
  expect_named(ps, c("arm1", "arm2", "average"))
  w <- colSums(f$p_z)
  for (nm in names(f$b)) {
    ts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    expect_equal(ps$average[[nm]],
                 (ps$arm1[[nm]] * w[[ts[1]]] + ps$arm2[[nm]] * w[[ts[2]]]) /
                   (w[[ts[1]]] + w[[ts[2]]]))
  }
})

test_that("solutions are named by the adherence sets they select", {
  f <- thesis_fit()
  # One solution per contrast here. t1 vs t3 is identified only for decision
  # teams whose adherence set is exactly {t1, t3}; t1 vs t2 for those whose
  # set contains both, that is A4 = {t1,t2} or A7 = {t1,t2,t3}.
  expect_identical(rownames(f$b$t1_t3), "A5")
  expect_identical(rownames(f$b$t1_t2), "A4+A7")
  ps <- estimate_p_sigma(f$p_z, f$KB, f$b)
  expect_identical(names(ps$arm1$t1_t3), "A5")
})

test_that("estimate_p_sigma recovers the true size of the target population", {
  f <- thesis_fit()
  ps <- estimate_p_sigma(f$p_z, f$KB, f$b)
  expect_equal(unname(ps$average$t1_t3["A5"]),
               mean(f$d$a == "A5"), tolerance = 0.02)
  expect_equal(unname(ps$average$t1_t2["A4+A7"]),
               mean(f$d$a %in% c("A4", "A7")), tolerance = 0.02)
})

test_that("estimate_late is the difference of the two arm responses", {
  f <- thesis_fit()
  ps <- estimate_p_sigma(f$p_z, f$KB, f$b)
  both <- estimate_late(f$q_z, f$KB, f$b, ps, arms = TRUE)
  expect_equal(both$effect, estimate_late(f$q_z, f$KB, f$b, ps))
  for (nm in names(both$effect)) {
    expect_equal(both$arm1[[nm]] - both$arm2[[nm]], both$effect[[nm]])
  }
})

test_that("the shared denominator cancels in a ratio", {
  f <- thesis_fit()
  ps <- estimate_p_sigma(f$p_z, f$KB, f$b)
  ratio_avg <- estimate_late(f$q_z, f$KB, f$b, ps,
                             scale = "ratio", denominator = "average")
  # With a common denominator the risk ratio is the ratio of the raw
  # projected quantities.
  for (nm in names(f$b)) {
    ts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    n1 <- as.numeric(f$b[[nm]] %*% (f$KB[[ts[1]]]$B_plus %*% f$q_z[, ts[1]]))
    n2 <- as.numeric(f$b[[nm]] %*% (f$KB[[ts[2]]]$B_plus %*% f$q_z[, ts[2]]))
    expect_equal(unname(ratio_avg[[nm]]), n1 / n2)
  }
})

test_that("estimated effects recover the simulated truth", {
  f <- thesis_fit()
  ps <- estimate_p_sigma(f$p_z, f$KB, f$b)
  late <- estimate_late(f$q_z, f$KB, f$b, ps)
  # Relative tolerance of 0.1. At n = 30000 the single-draw Monte Carlo error
  # on these ratio estimators is around 0.01 in absolute terms, i.e. about 5%
  # of a 0.2 contrast. Averaged over seeds the estimates converge on the true
  # values as n grows, so this margin reflects sampling noise, not bias.
  expect_equal(unname(late$t1_t2["A4+A7"]),
               unname(f$truth["t1"] - f$truth["t2"]), tolerance = 0.1)
  expect_equal(unname(late$t2_t3["A6+A7"]),
               unname(f$truth["t2"] - f$truth["t3"]), tolerance = 0.1)
  expect_equal(unname(late$t1_t3["A5"]),
               unname(f$truth["t1"] - f$truth["t3"]), tolerance = 0.1)
})

test_that("the IV estimate beats the confounded comparison", {
  f <- thesis_fit()
  late <- estimate_late(f$q_z, f$KB, f$b,
                        estimate_p_sigma(f$p_z, f$KB, f$b))
  target <- unname(f$truth["t1"] - f$truth["t3"])
  crude <- mean(f$d$y[f$d$trt == "t1"]) - mean(f$d$y[f$d$trt == "t3"])
  expect_gt(abs(crude - target), abs(unname(late$t1_t3["A5"]) - target))
})

test_that("the single-treatment estimators give responses, not contrasts", {
  f <- thesis_fit()
  b_t <- solve_b_treatments(f$KB)
  ps <- estimate_p_sigma_treatments(f$p_z, f$KB, b_t)
  latr <- estimate_latr(f$q_z, f$KB, b_t, ps)
  expect_named(latr, f$treatments)
  # Responses are contaminated by the adherence-set shift, but their
  # differences are not: contrasts of responses reproduce the LATE.
  late <- estimate_late(f$q_z, f$KB, f$b,
                        estimate_p_sigma(f$p_z, f$KB, f$b))
  expect_equal(unname(latr$t1["A4+A7"] - latr$t2["A4+A7"]),
               unname(late$t1_t2["A4+A7"]))
  expect_equal(unname(latr$t1["A5"] - latr$t3["A5"]),
               unname(late$t1_t3["A5"]))
})
