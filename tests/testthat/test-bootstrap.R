boot_fit <- function(n = 1500, seed = 3) {
  treatments <- c("t1", "t2", "t3")
  instrument <- list(`1` = c("t1", "t2", "t3"), `2` = c("t2", "t1", "t3"),
                     `3` = c("t2", "t3", "t1"), `4` = c("t3", "t2", "t1"))
  R <- response_matrix(adherence_sets(treatments), instrument)
  KB <- projection_matrices(R, treatments)
  sets <- adherence_sets(treatments)
  set.seed(seed)
  a <- sample(names(sets), n, TRUE, prob = c(1, 1, 1, 2, 3, 2, 4))
  z <- sample(seq_along(instrument), n, TRUE)
  trt <- vapply(seq_len(n),
                function(i) choose_treatment(instrument[[z[i]]], sets[[a[i]]]),
                character(1))
  d <- data.frame(z = z, trt = trt, x = stats::rnorm(n),
                  y = stats::rbinom(n, 1, c(t1 = .3, t2 = .5, t3 = .7)[trt]))
  p_z <- estimate_p_z(d, "z", "trt", instrument_levels = rownames(R),
                      treatments = treatments)
  list(treatments = treatments, instrument = instrument, R = R, KB = KB, d = d,
       b = solve_b_pairs(KB), p_z = p_z,
       q_z = estimate_q_z(d, "z", "trt", "y", instrument_levels = rownames(R),
                          treatments = treatments),
       p_sigma = estimate_p_sigma(p_z, KB, solve_b_pairs(KB)))
}

test_that("a weighted comparison on the pseudo-population reproduces the LATE", {
  f <- boot_fit()
  late <- estimate_late(f$q_z, f$KB, f$b, f$p_sigma)
  for (nm in names(f$b)) {
    ts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    for (sol in rownames(f$b[[nm]])) {
      p <- pseudo_population(nm, f$d, f$KB, f$b, f$p_sigma, "z", "trt",
                             solution = sol)
      arm1 <- as.character(p$trt) == ts[1]
      expect_equal(
        weighted.mean(p$y[arm1], p$w[arm1]) -
          weighted.mean(p$y[!arm1], p$w[!arm1]),
        unname(late[[nm]][[sol]])
      )
    }
  }
})

test_that("pseudo-population weights average to one within each arm", {
  # This is what makes mean(x * w) equal the weighted mean of x, which is how
  # the baseline tables in paper 2 are built.
  f <- boot_fit()
  for (nm in names(f$b)) {
    ts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    p <- pseudo_population(nm, f$d, f$KB, f$b, f$p_sigma, "z", "trt")
    for (t in ts) {
      arm <- as.character(p$trt) == t
      expect_equal(mean(p$w[arm]), 1)
      expect_equal(mean(p$x[arm] * p$w[arm]), weighted.mean(p$x[arm], p$w[arm]))
    }
  }
})

test_that("pseudo_population keeps only the two treatments compared", {
  f <- boot_fit()
  p <- pseudo_population("t1_t2", f$d, f$KB, f$b, f$p_sigma, "z", "trt")
  expect_setequal(unique(as.character(p$trt)), c("t1", "t2"))
  expect_true("w" %in% names(p))
  expect_equal(nrow(p), sum(f$d$trt %in% c("t1", "t2")))
})

test_that("solutions can be selected by name or by index, and bad ones error", {
  f <- boot_fit()
  by_index <- pseudo_population("t1_t2", f$d, f$KB, f$b, f$p_sigma, "z", "trt",
                                solution = 1L)
  by_name <- pseudo_population("t1_t2", f$d, f$KB, f$b, f$p_sigma, "z", "trt",
                               solution = rownames(f$b$t1_t2)[1])
  expect_equal(by_index$w, by_name$w)
  expect_error(
    pseudo_population("t1_t2", f$d, f$KB, f$b, f$p_sigma, "z", "trt",
                      solution = "A1"),
    "not identified"
  )
  expect_error(
    pseudo_population("t1_t2", f$d, f$KB, f$b, f$p_sigma, "z", "trt",
                      solution = 99L),
    "must be between"
  )
})

test_that("naive_iv is the Wald ratio on the binarised instrument", {
  f <- boot_fit()
  got <- naive_iv("t1_t2", f$d, f$instrument, "z", "trt", "y")
  # Only instrument value 1, c("t1", "t2", "t3"), puts t1 ahead of t2.
  ahead <- which(vapply(f$instrument,
                        function(z) which(z == "t1") < which(z == "t2"),
                        logical(1)))
  expect_equal(unname(ahead), 1L)
  d <- f$d[f$d$trt %in% c("t1", "t2"), ]
  first <- d$z %in% ahead
  expect_equal(
    got,
    (mean(d$y[first]) - mean(d$y[!first])) /
      (mean(d$trt[first] == "t1") - mean(d$trt[!first] == "t1"))
  )
})

test_that("bootstrap_late returns intervals that bracket the point estimate", {
  skip_on_cran()
  f <- boot_fit()
  bs <- bootstrap_late(24, f$d, "z", "trt", "y", f$KB, f$b, n_cores = 2)
  late <- estimate_late(f$q_z, f$KB, f$b, f$p_sigma)
  expect_named(bs, c("alpha", "replicates", "ci"))
  expect_equal(nrow(bs$replicates), 24L)
  expect_true(all(bs$ci$lower <= bs$ci$upper))
  for (nm in names(late)) {
    for (sol in names(late[[nm]])) {
      key <- paste(nm, sol, sep = ".")
      expect_true(key %in% rownames(bs$ci))
      expect_gte(late[[nm]][[sol]], bs$ci[key, "lower"] - 1e-8)
      expect_lte(late[[nm]][[sol]], bs$ci[key, "upper"] + 1e-8)
    }
  }
})

test_that("bootstrap_late is reproducible given the same seed and worker count", {
  skip_on_cran()
  f <- boot_fit()
  a <- bootstrap_late(16, f$d, "z", "trt", "y", f$KB, f$b, n_cores = 2)
  b2 <- bootstrap_late(16, f$d, "z", "trt", "y", f$KB, f$b, n_cores = 2)
  expect_equal(a$replicates, b2$replicates)
  c3 <- bootstrap_late(16, f$d, "z", "trt", "y", f$KB, f$b, n_cores = 2,
                       seed = 99)
  expect_false(isTRUE(all.equal(a$replicates, c3$replicates)))
})

test_that("a single identified effect keeps its shape", {
  skip_on_cran()
  # With two treatments there is exactly one contrast and one solution, which
  # is the case where a collapsing sapply would lose the row per effect.
  treatments <- c("a", "b")
  instrument <- list(`1` = c("a", "b"), `2` = c("b", "a"))
  R <- response_matrix(adherence_sets(treatments), instrument)
  KB <- projection_matrices(R, treatments)
  b <- solve_b_pairs(KB)
  set.seed(8)
  d <- data.frame(
    z = rep(1:2, each = 200),
    trt = c(sample(c("a", "b"), 200, TRUE, c(0.7, 0.3)),
            sample(c("a", "b"), 200, TRUE, c(0.3, 0.7))),
    y = stats::rbinom(400, 1, 0.5)
  )
  bs <- bootstrap_late(12, d, "z", "trt", "y", KB, b, n_cores = 2)
  expect_equal(dim(bs$replicates), c(12L, 1L))
  expect_equal(nrow(bs$ci), 1L)
  expect_identical(rownames(bs$ci), "a_b.A3")
})

test_that("a list of imputations gives n replicates per imputation", {
  skip_on_cran()
  f <- boot_fit()
  imps <- list(f$d, f$d, f$d)
  bs <- bootstrap_late(8, imps, "z", "trt", "y", f$KB, f$b, n_cores = 2)
  expect_equal(nrow(bs$replicates), 24L)
  expect_error(bootstrap_late(8, list(), "z", "trt", "y", f$KB, f$b,
                              n_cores = 2),
               "data frame")
})
