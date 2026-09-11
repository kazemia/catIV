# The worked example from the thesis, Summary of Simulations section 3.1:
# three treatments, four instrument values. Three effects are identifiable,
# with the target populations stated there.
thesis_setup <- function() {
  treatments <- c("t1", "t2", "t3")
  instrument <- list(z1 = c("t1", "t2", "t3"), z2 = c("t2", "t1", "t3"),
                     z3 = c("t2", "t3", "t1"), z4 = c("t3", "t2", "t1"))
  R <- response_matrix(adherence_sets(treatments), instrument)
  list(treatments = treatments, instrument = instrument, R = R,
       KB = projection_matrices(R, treatments))
}

test_that("adherence_sets enumerates 2^NT - 1 non-empty sets", {
  expect_length(adherence_sets(c("a", "b", "c")), 7L)
  expect_length(adherence_sets(letters[1:4]), 15L)
  expect_identical(adherence_sets(c("a", "b"))$A1, "a")
  expect_error(adherence_sets(c("a", "a")), "duplicates")
  expect_error(adherence_sets("a"), "at least two")
})

test_that("choose_treatment picks the most encouraged adherent alternative", {
  expect_identical(choose_treatment(c("a", "b", "c"), c("b", "c")), "b")
  expect_identical(choose_treatment(c("c", "b", "a"), c("b", "c")), "c")
  expect_identical(choose_treatment(c("a", "b", "c"), "c"), "c")
  expect_error(choose_treatment(c("a", "b"), "z"), "No element")
})

test_that("response_matrix has one row per instrument value and one column per set", {
  s <- thesis_setup()
  expect_equal(dim(s$R), c(4L, 7L))
  expect_identical(rownames(s$R), names(s$instrument))
  expect_identical(colnames(s$R), names(adherence_sets(s$treatments)))
  # A singleton adherence set always receives its own treatment.
  expect_true(all(s$R[["A1"]] == "t1"))
})

test_that("K is the projector onto the null space of B", {
  s <- thesis_setup()
  for (t in s$treatments) {
    B <- s$KB[[t]]$B
    K <- s$KB[[t]]$K
    # b K = 0 is equivalent to b = b B^+ B, so K must annihilate the row space.
    expect_equal(unname(B %*% K), matrix(0, nrow(B), ncol(K)), tolerance = 1e-6)
    expect_equal(unname(K %*% K), unname(K), tolerance = 1e-6)   # idempotent
    expect_equal(unname(t(K)), unname(K), tolerance = 1e-6)      # symmetric
  }
})

test_that("solve_b_pairs reproduces the identifiable effects in the thesis", {
  s <- thesis_setup()
  pops <- target_populations(solve_b_pairs(s$KB))
  # Adherence sets: A1={t1} A2={t2} A3={t3} A4={t1,t2} A5={t1,t3}
  #                 A6={t2,t3} A7={t1,t2,t3}
  expect_setequal(unlist(pops$t1_t2), c("A4", "A7"))  # {t1,t2} subset of A
  expect_setequal(unlist(pops$t2_t3), c("A6", "A7"))  # {t2,t3} subset of A
  expect_setequal(unlist(pops$t1_t3), "A5")           # A exactly {t1,t3}
})

test_that("every returned b actually solves the defining equation", {
  s <- thesis_setup()
  b_pairs <- solve_b_pairs(s$KB)
  for (nm in names(b_pairs)) {
    ts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    for (r in seq_len(nrow(b_pairs[[nm]]))) {
      bb <- b_pairs[[nm]][r, ]
      expect_lt(max(abs(bb %*% s$KB[[ts[1]]]$K)), 1e-3)
      expect_lt(max(abs(bb %*% s$KB[[ts[2]]]$K)), 1e-3)
      expect_gt(sum(bb), 0)
    }
  }
})

test_that("solve_b_treatments finds every binary solution, verified by brute force", {
  s <- thesis_setup()
  b_trt <- solve_b_treatments(s$KB)
  n <- ncol(s$R)
  grid <- as.matrix(expand.grid(rep(list(0:1), n)))
  colnames(grid) <- colnames(s$R)
  for (t in s$treatments) {
    K <- s$KB[[t]]$K
    ok <- apply(grid, 1, function(bb) sum(bb) > 0 && max(abs(bb %*% K)) < 1e-3)
    truth <- grid[ok, , drop = FALSE]
    got <- b_trt[[t]]
    key <- function(m) sort(unname(apply(m, 1, paste, collapse = "")))
    expect_identical(key(got), key(truth))
  }
})

test_that("solve_b_pairs is the simultaneous solution set of both treatments", {
  s <- thesis_setup()
  b_trt <- solve_b_treatments(s$KB)
  b_pairs <- solve_b_pairs(s$KB)
  key <- function(m) apply(m, 1, paste, collapse = "")
  for (nm in names(b_pairs)) {
    ts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    expect_setequal(key(b_pairs[[nm]]),
                    intersect(key(b_trt[[ts[1]]]), key(b_trt[[ts[2]]])))
  }
})

test_that("a pair with no identifiable effect returns a zero-row matrix", {
  # A single instrument value cannot identify anything.
  treatments <- c("a", "b")
  R <- response_matrix(adherence_sets(treatments), list(z1 = c("a", "b")))
  b <- solve_b_pairs(projection_matrices(R, treatments))
  expect_equal(nrow(b$a_b), 0L)
  expect_identical(colnames(b$a_b), c("A1", "A2", "A3"))
})
