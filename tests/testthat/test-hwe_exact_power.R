
test_that("Power calculation for theta = 8 is correct", {
  N <- 100
  K <- 14  # Equivalent to nA
  theta <- 8
  expected_power <- 0.0008
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE)
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})

test_that("Power calculation for theta = 4 is correct", {
  N <- 100
  K <- 14
  theta <- 4
  expected_power <- 0.0053
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE)
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})

test_that("Power calculation for theta = 2 is correct", {
  N <- 100
  K <- 14
  theta <- 2
  expected_power <- 0.0285
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE)
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})

test_that("Power calculation for theta = 0.1 is correct", {
  N <- 100
  K <- 14
  theta <- 0.1
  expected_power <- 0.9185
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE)
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})


test_that("Power calculation for theta = 8 is correct when alternative is 'greater'", {
  N <- 100
  K <- 54  # Equivalent to nA
  theta <- 8
  expected_power <- 0.247597
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE, alternative = "greater")
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})

test_that("Power calculation for theta = 4 is correct when alternative is 'greater'", {
  N <- 100
  K <- 54
  theta <- 4
  expected_power <- 0.0254652
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE, alternative = "greater")
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})

test_that("Power calculation for theta = 2 is correct when alternative is 'greater'", {
  N <- 100
  K <- 54
  theta <- 2
  expected_power <- 0.0004663818
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE, alternative = "greater")
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})

test_that("Power calculation for theta = 0.1 is correct when alternative is 'greater'", {
  N <- 100
  K <- 54
  theta <- 0.1
  expected_power <- 4.934875e-21
  calculated_power <- hwe_exact_power(N, K, theta, sig.level = 0.05, use_mid_p = FALSE, alternative = "greater")
  expect_equal(calculated_power, expected_power, tolerance = 1e-4)
})



