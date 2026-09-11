test_that("the three principal strata partition the population", {
  tr <- binary_iv_truth(seq(-3, 3, by = 0.5), z_on_t = 1.5, x_on_t = 0.4,
                        u_on_t = 1, v_on_t = 0.8, u_prob = 0.5, v_prob = 0.6)
  expect_equal(tr$complier + tr$always_taker + tr$never_taker,
               rep(1, nrow(tr)))
  expect_true(all(tr$complier >= 0))
  expect_true(all(tr$always_taker >= 0 & tr$always_taker <= 1))
  expect_true(all(tr$never_taker >= 0 & tr$never_taker <= 1))
})

test_that("a stronger instrument makes more compliers", {
  weak <- binary_iv_truth(0, z_on_t = 0.2, x_on_t = 0, u_on_t = 0, v_on_t = 0,
                          u_prob = 0.5, v_prob = 0.5)
  strong <- binary_iv_truth(0, z_on_t = 3, x_on_t = 0, u_on_t = 0, v_on_t = 0,
                            u_prob = 0.5, v_prob = 0.5)
  expect_gt(strong$complier, weak$complier)
  # With no instrument at all there are no compliers.
  none <- binary_iv_truth(0, z_on_t = 0, x_on_t = 0, u_on_t = 0, v_on_t = 0,
                          u_prob = 0.5, v_prob = 0.5)
  expect_equal(none$complier, 0)
})

test_that("the marginal treated probability lies between the two arms", {
  tr <- binary_iv_truth(c(-1, 0, 1), z_on_t = 1.5, x_on_t = 0.4, u_on_t = 1,
                        v_on_t = 0.8, u_prob = 0.5, v_prob = 0.6, z_prob = 0.3)
  encouraged <- 1 - tr$never_taker
  expect_true(all(tr$treated >= tr$always_taker & tr$treated <= encouraged))
  # z_prob shifts it between the two arms linearly.
  expect_equal(tr$treated, 0.3 * encouraged + 0.7 * tr$always_taker)
})

test_that("the observed confounder is marginalised over", {
  args <- list(x = c(0, 1), z_on_t = 1.5, x_on_t = 0.4, u_on_t = 1,
               u_prob = 0.5, v_prob = 0.6)
  with_v <- do.call(binary_iv_truth, c(args, list(v_on_t = 0.8)))
  without_v <- do.call(binary_iv_truth, c(args, list(v_on_t = 0)))
  # Every quantity, the marginal one included, must respond to the V effect.
  expect_false(isTRUE(all.equal(with_v$treated, without_v$treated)))
  expect_false(isTRUE(all.equal(with_v$always_taker, without_v$always_taker)))
  # With v_prob = 0 the confounder never occurs, so its effect cannot matter.
  a <- binary_iv_truth(0, 1.5, 0.4, 1, v_on_t = 0.8, u_prob = 0.5, v_prob = 0)
  b <- binary_iv_truth(0, 1.5, 0.4, 1, v_on_t = 0, u_prob = 0.5, v_prob = 0)
  expect_equal(a, b)
})

test_that("results are vectorised over the covariate", {
  xs <- seq(-2, 2, length.out = 7)
  tr <- binary_iv_truth(xs, 1.5, 0.4, 1, 0.8, 0.5, 0.5)
  expect_equal(nrow(tr), length(xs))
  expect_equal(tr$x, xs)
  one <- binary_iv_truth(xs[3], 1.5, 0.4, 1, 0.8, 0.5, 0.5)
  expect_equal(as.numeric(tr[3, ]), as.numeric(one[1, ]))
})
