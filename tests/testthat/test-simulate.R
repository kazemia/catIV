test_that("simulate_binary returns binary columns of the right length", {
  set.seed(1)
  d <- simulate_binary(200, z_on_t = 0.4, v_on_t = 0.3, v_on_y = 0.25,
                       t_on_y = 0.2)
  expect_named(d, c("Z", "T", "Y"))
  expect_equal(nrow(d), 200L)
  expect_true(all(vapply(d, function(x) all(x %in% c(0, 1)), logical(1))))
})

test_that("simulate_binary respects the instrument and treatment effects", {
  set.seed(2)
  d <- simulate_binary(20000, z_on_t = 0.4, v_on_t = 0.3, v_on_y = 0.25,
                       t_on_y = 0.2)
  # The instrument shifts treatment uptake by z_on_t.
  expect_equal(mean(d$T[d$Z == 1]) - mean(d$T[d$Z == 0]), 0.4, tolerance = 0.05)
  # Treatment raises the outcome, though confounded upwards by V.
  expect_gt(mean(d$Y[d$T == 1]) - mean(d$Y[d$T == 0]), 0.2)
})

test_that("the logit link keeps probabilities valid for large effects", {
  set.seed(3)
  d <- simulate_binary(500, z_on_t = 4, v_on_t = 3, v_on_y = 2, t_on_y = 5,
                       logit = TRUE)
  expect_false(anyNA(d))
  # On the additive scale the same effects put probabilities outside [0, 1],
  # which rbinom turns into NA. This is the caller's responsibility, and is
  # why `logit` exists.
  bad <- suppressWarnings(
    simulate_binary(500, z_on_t = 4, v_on_t = 3, v_on_y = 2, t_on_y = 5,
                    logit = FALSE)
  )
  expect_true(anyNA(bad))
})

sim_setup <- function() {
  treatments <- c("t1", "t2", "t3")
  instrument <- list(`1` = c("t1", "t2", "t3"), `2` = c("t2", "t1", "t3"),
                     `3` = c("t2", "t3", "t1"), `4` = c("t3", "t2", "t1"))
  sets <- adherence_sets(treatments)
  set.seed(9)
  v_to_set <- lapply(list("0" = NULL, "1" = NULL),
                     function(z) { p <- stats::runif(length(sets)); p / sum(p) })
  list(treatments = treatments, instrument = instrument, sets = sets,
       v_to_set = v_to_set, v_on_y = c(0.10, 0.25),
       t_on_y = c(0.20, 0.35, 0.50))
}

test_that("simulate_adherence obeys categorical monotonicity exactly", {
  s <- sim_setup()
  set.seed(11)
  d <- simulate_adherence(1000, s$instrument, s$sets, s$treatments,
                          v_levels = c(0, 1), v_probs = c(0.5, 0.5),
                          v_to_set = s$v_to_set, v_on_y = s$v_on_y,
                          t_on_y = s$t_on_y)
  expect_named(d, c("Z", "T", "Y", "A"))
  # Every treatment received is the one the choice function picks.
  implied <- vapply(seq_len(nrow(d)),
                    function(i) choose_treatment(s$instrument[[d$Z[i]]],
                                                 s$sets[[d$A[i]]]),
                    character(1))
  expect_identical(d$T, implied)
  # The treatment received always belongs to the adherence set.
  expect_true(all(vapply(seq_len(nrow(d)),
                         function(i) d$T[i] %in% s$sets[[d$A[i]]],
                         logical(1))))
})

test_that("a singleton adherence set always takes its own treatment", {
  s <- sim_setup()
  only_t2 <- rep(0, length(s$sets))
  only_t2[match("A2", names(s$sets))] <- 1   # A2 = {t2}
  set.seed(12)
  d <- simulate_adherence(200, s$instrument, s$sets, s$treatments,
                          v_levels = 0, v_probs = 1,
                          v_to_set = list("0" = only_t2),
                          v_on_y = 0, t_on_y = s$t_on_y)
  expect_true(all(d$A == "A2"))
  expect_true(all(d$T == "t2"))
})

test_that("simulate_categorical returns only the observed confounders", {
  treatments <- c("t1", "t2", "t3")
  instrument <- list(`1` = c("t1", "t2", "t3"), `2` = c("t3", "t2", "t1"))
  v_on_t <- matrix(c(3, 0, 0, 3, 1, 1), nrow = 3, byrow = TRUE)
  args <- list(n = 400, instrument = instrument, treatments = treatments,
               v_on_t = v_on_t, v_probs = c(0.4, 0.6), v_on_y = c(0.12, 0.08),
               t_on_y = c(0.25, 0.35, 0.45))
  set.seed(4); both <- do.call(simulate_categorical, args)
  set.seed(4); one <- do.call(simulate_categorical,
                              c(args, list(observed = c(TRUE, FALSE))))
  expect_named(both, c("Z", "T", "Y", "V1", "V2"))
  expect_named(one, c("Z", "T", "Y", "V1"))
  # Hiding a confounder must not change the data that is returned.
  expect_identical(both[c("Z", "T", "Y", "V1")], one)
})

test_that("simulate_categorical keeps outcome probabilities valid", {
  treatments <- c("t1", "t2")
  instrument <- list(`1` = c("t1", "t2"), `2` = c("t2", "t1"))
  set.seed(6)
  # Effects that would exceed 1 are clipped rather than erroring.
  d <- simulate_categorical(300, instrument, treatments,
                            v_on_t = matrix(c(2, 0), nrow = 2),
                            v_probs = 0.5, v_on_y = 0.9, t_on_y = c(0.5, 0.8))
  expect_true(all(d$Y %in% c(0, 1)))
  expect_setequal(unique(d$T), treatments)
})

test_that("the instrument pushes treatment in the encouraged direction", {
  treatments <- c("t1", "t2", "t3")
  instrument <- list(`1` = c("t1", "t2", "t3"), `2` = c("t3", "t2", "t1"))
  set.seed(7)
  d <- simulate_categorical(8000, instrument, treatments,
                            v_on_t = matrix(0, nrow = 3, ncol = 1),
                            v_probs = 0.5, v_on_y = 0, t_on_y = c(.2, .3, .4))
  expect_gt(mean(d$T[d$Z == 1] == "t1"), mean(d$T[d$Z == 2] == "t1"))
  expect_gt(mean(d$T[d$Z == 2] == "t3"), mean(d$T[d$Z == 1] == "t3"))
})
