test_that("multivariate independent components works correctly", {
  mu1 <- c(0, 1, 2)
  sigma1 <- c(1, 1.5, 2)
  nu1 <- c(10, 12, 15)
  mu2 <- c(0, 0, 1)
  sigma2 <- c(1.2, 1, 1.8)
  nu2 <- c(15, 20, 25)

  result <- ws_tdiff_multivariate_independent(mu1, sigma1, nu1, mu2, sigma2, nu2)

  # Check structure
  expect_s3_class(result, "ws_tdiff_multivariate_independent")
  expect_equal(length(result$mu_diff), 3)
  expect_equal(length(result$sigma_star), 3)
  expect_equal(length(result$nu_star), 3)
  expect_equal(result$p, 3)

  # Check component-wise matching
  for (j in 1:3) {
    uni_result <- ws_tdiff_univariate(
      mu1[j], sigma1[j], nu1[j],
      mu2[j], sigma2[j], nu2[j]
    )
    expect_equal(result$mu_diff[j], uni_result$mu_diff)
    expect_equal(result$sigma_star[j], uni_result$sigma_star)
    expect_equal(result$nu_star[j], uni_result$nu_star)
  }
})

test_that("multivariate independent validates dimensions", {
  expect_error(
    ws_tdiff_multivariate_independent(
      c(0, 1), c(1), c(10, 12),  # Mismatched dimensions
      c(0, 0), c(1.2, 1), c(15, 20)
    ),
    "same length"
  )

  expect_error(
    ws_tdiff_multivariate_independent(
      c(0, 1), c(1, 1.5), c(3, 12),  # nu1[1] <= 4
      c(0, 0), c(1.2, 1), c(15, 20)
    ),
    "greater than 4"
  )
})

test_that("multivariate general case works correctly", {
  # Create positive definite matrices
  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)

  result <- ws_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  # Check structure
  expect_s3_class(result, "ws_tdiff_multivariate_general")
  expect_true(is.matrix(result$Sigma_star))
  expect_equal(dim(result$Sigma_star), c(2, 2))
  expect_true(result$nu_star > 0)
  expect_true(result$converged)

  # Check symmetry of result
  expect_true(isSymmetric(result$Sigma_star))

  # Check positive definiteness
  eigenvalues <- eigen(result$Sigma_star, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("multivariate general validates input", {
  # Not a matrix
  expect_error(
    ws_tdiff_multivariate_general(
      c(0, 1), c(1, 0, 0, 1), 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "must be matrices"
  )

  # Non-square matrix
  expect_error(
    ws_tdiff_multivariate_general(
      c(0, 1), matrix(1:6, 2, 3), 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "square"
  )

  # Dimension mismatch
  expect_error(
    ws_tdiff_multivariate_general(
      c(0, 1, 2), matrix(c(1, 0, 0, 1), 2), 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "consistent"
  )

  # Non-positive definite
  bad_matrix <- matrix(c(1, 2, 2, 1), 2, 2)  # Not positive definite
  expect_error(
    ws_tdiff_multivariate_general(
      c(0, 1), bad_matrix, 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "positive definite"
  )
})

test_that("diagonal matrices reduce to independent case", {
  # When using diagonal matrices, should get similar results
  Sigma1 <- diag(c(1, 1.5))
  Sigma2 <- diag(c(1.2, 1))

  result_general <- ws_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  result_indep <- ws_tdiff_multivariate_independent(
    mu1 = c(0, 1), sigma1 = c(1, 1.5), nu1 = c(10, 10),
    mu2 = c(0, 0), sigma2 = c(1.2, 1), nu2 = c(15, 15)
  )

  # Location should match exactly
  expect_equal(result_general$mu_diff, result_indep$mu_diff)

  # Diagonal elements should be close
  diag_general <- diag(result_general$Sigma_star)
  # Note: These won't match exactly due to different df treatment
  expect_true(all(abs(sqrt(diag_general) - result_indep$sigma_star) < 1))
})
