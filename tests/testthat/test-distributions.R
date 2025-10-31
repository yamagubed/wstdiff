test_that("distribution functions work correctly", {
  result <- ws_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # Test density
  d <- dtdiff(0, result)
  expect_true(is.numeric(d))
  expect_true(d > 0)
  expect_true(is.finite(d))

  # Test CDF
  p <- ptdiff(0, result)
  expect_true(p >= 0 && p <= 1)
  expect_equal(ptdiff(-Inf, result), 0)
  expect_equal(ptdiff(Inf, result), 1)

  # Test monotonicity
  x_seq <- seq(-5, 5, by = 0.5)
  p_seq <- ptdiff(x_seq, result)
  expect_true(all(diff(p_seq) >= 0))  # Should be non-decreasing

  # Test quantile
  q <- qtdiff(0.5, result)
  expect_true(is.numeric(q))
  expect_equal(q, result$mu_diff)  # Median should be mu_diff

  # Test quantile boundaries
  expect_equal(qtdiff(0, result), -Inf)
  expect_equal(qtdiff(1, result), Inf)

  # Test random generation
  set.seed(123)
  r <- rtdiff(1000, result)
  expect_equal(length(r), 1000)
  expect_true(all(is.finite(r)))

  # Mean should be close to mu_diff
  expect_equal(mean(r), result$mu_diff, tolerance = 0.1)
})

test_that("distribution functions validate input", {
  result <- ws_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # Invalid probability for qtdiff
  expect_error(qtdiff(1.5, result), "between 0 and 1")
  expect_error(qtdiff(-0.1, result), "between 0 and 1")

  # Invalid ws_result
  expect_error(dtdiff(0, list()), "ws_tdiff_univariate")
  expect_error(ptdiff(0, "not_a_result"), "ws_tdiff_univariate")
})

test_that("ptdiff and qtdiff are inverses", {
  result <- ws_tdiff_univariate(2, 1.5, 12, 1, 2, 18)

  # Test multiple probability points
  probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

  for (p in probs) {
    q <- qtdiff(p, result)
    p_back <- ptdiff(q, result)
    expect_equal(p, p_back, tolerance = 1e-10)
  }

  # Test multiple quantile points
  quantiles <- seq(-5, 5, by = 1)
  for (q in quantiles) {
    p <- ptdiff(q, result)
    q_back <- qtdiff(p, result)
    expect_equal(q, q_back, tolerance = 1e-10)
  }
})

test_that("density integrates to 1", {
  result <- ws_tdiff_univariate(0, 1, 10, 0, 1, 10)

  # Numerical integration of density
  integrate_density <- integrate(
    function(x) dtdiff(x, result),
    lower = -Inf,
    upper = Inf
  )

  expect_equal(integrate_density$value, 1, tolerance = 1e-6)
})
