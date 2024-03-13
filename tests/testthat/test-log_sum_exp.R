
test_that("Log-Sum-Exp returns correct value with finite numbers", {
  log_values <- c(log(0.1), log(0.2), log(0.7))
  expect_equal(log_sum_exp(log_values), log(1), tolerance = 1e-8)
})

test_that("Log-Sum-Exp handles -Inf values", {
  log_values <- c(-Inf, log(0.1), log(0.9))
  expect_equal(log_sum_exp(log_values), log(1), tolerance = 1e-8)
})

test_that("Log-Sum-Exp returns -Inf with all -Inf values", {
  log_values <- rep(-Inf, 3)
  expect_equal(log_sum_exp(log_values), -Inf)
})

test_that("Log-Sum-Exp handles log(0) correctly", {
  log_values <- c(log(0), log(0), log(0.5))
  expect_equal(log_sum_exp(log_values), log(0.5), tolerance = 1e-8)
})


test_that("log_sum_exp correctly sums log-transformed values", {
  log_values <- log(c(1, 10, 100))
  expect_equal(log_sum_exp(log_values), log(sum(exp(log_values))))
})

test_that("log_sum_exp handles -Inf correctly", {
  log_values <- c(log(10), -Inf, log(100))
  expect_equal(log_sum_exp(log_values), log(sum(exp(log_values))))
})

test_that("log_sum_exp is accurate with large differences", {
  log_values <- c(log(1e-100), log(1e+100))
  expect_equal(log_sum_exp(log_values), log(sum(exp(log_values))))
})

test_that("log_sum_exp returns -Inf for empty input", {
  log_values <- numeric(0) # Empty vector
  expect_equal(log_sum_exp(log_values), -Inf)
})

test_that("log_sum_exp returns the same value for single input", {
  log_value <- log(100)
  expect_equal(log_sum_exp(c(log_value)), log_value)
})

test_that("log_sum_exp correctly computes the log of sums of exponentials", {
  # Test with a straightforward case
  expect_equal(log_sum_exp(c(log(2), log(3))), log(5))
  
  # Test with large numbers that could cause overflow in a naive implementation
  large_values <- c(log(1e300), log(2e300))
  expect_equal(log_sum_exp(large_values), log(3e300))
  
  # Test with small numbers that could cause underflow in a naive implementation
  small_values <- c(log(1e-300), log(2e-300))
  expect_equal(log_sum_exp(small_values), log(3e-300))
  
  # Test with a mix of large, small, and -Inf values
  mixed_values <- c(log(1e300), log(1e-300), -Inf)
  expect_equal(log_sum_exp(mixed_values), log(1e300 + 1e-300))
  
  # Test with only -Inf values (expect -Inf as result)
  inf_values <- rep(-Inf, 5)
  expect_equal(log_sum_exp(inf_values), -Inf)
  
  # Test with an empty vector (expect -Inf as result)
  expect_equal(log_sum_exp(numeric(0)), -Inf)
})


