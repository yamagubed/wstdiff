test_that("validation function works", {
  result <- ws_tdiff_univariate(0, 1, 10, 0, 1.5, 15)
  validation <- validate_approximation(result)

  expect_type(validation, "list")
  expect_true("mean" %in% names(validation))
  expect_true("variance" %in% names(validation))
  expect_true("effective_df" %in% names(validation))

  # Check computed values
  expect_equal(validation$mean, result$mu_diff)
  expect_equal(validation$effective_df, result$nu_star)

  # Variance should match
  theoretical_var <- result$sigma_star^2 * result$nu_star / (result$nu_star - 2)
  expect_equal(validation$variance, theoretical_var)
})

test_that("validation provides appropriate warnings", {
  # Low degrees of freedom
  result_low <- ws_tdiff_univariate(0, 1, 5, 0, 1, 6)
  validation_low <- validate_approximation(result_low)
  expect_true(grepl("Warning", validation_low$quality_note))

  # High degrees of freedom
  result_high <- ws_tdiff_univariate(0, 1, 50, 0, 1, 60)
  validation_high <- validate_approximation(result_high)
  expect_true(grepl("good", validation_high$quality_note))
})
