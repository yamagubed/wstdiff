test_that("equal parameters formula is mathematically correct", {
  # Test various parameter combinations
  test_cases <- list(
    list(mu = 0, sigma = 1, nu = 10),
    list(mu = 5, sigma = 2, nu = 20),
    list(mu = -3, sigma = 0.5, nu = 100)
  )

  for (case in test_cases) {
    result <- ws_tdiff_equal_params(case$mu, case$sigma, case$nu)

    # Check formula: nu_star = 2(nu - 4)
    expect_equal(result$nu_star, 2 * (case$nu - 4))

    # Check formula: sigma_star = sigma * sqrt(2*nu/(nu-2))
    expected_sigma <- case$sigma * sqrt(2 * case$nu / (case$nu - 2))
    expect_equal(result$sigma_star, expected_sigma)

    # Mean difference should be 0
    expect_equal(result$mu_diff, 0)
  }
})

test_that("equal parameters matches general formula", {
  # Multiple test cases
  params <- list(
    list(mu = 0, sigma = 1, nu = 10),
    list(mu = 2, sigma = 3, nu = 50),
    list(mu = -1, sigma = 0.5, nu = 8)
  )

  for (p in params) {
    general <- ws_tdiff_univariate(
      p$mu, p$sigma, p$nu,
      p$mu, p$sigma, p$nu
    )
    special <- ws_tdiff_equal_params(p$mu, p$sigma, p$nu)

    expect_equal(general$mu_diff, special$mu_diff, tolerance = 1e-10)
    expect_equal(general$sigma_star, special$sigma_star, tolerance = 1e-10)
    expect_equal(general$nu_star, special$nu_star, tolerance = 1e-10)
  }
})
