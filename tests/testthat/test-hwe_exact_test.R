

test_that("hwe_exact_test calculates correct p-values", {
  # Test case 1
  N <- 100
  nA <- 21
  nAB <- 19
  theta <- 4
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "less", use_mid_p = TRUE, theta = theta),
    0.4872189,
    tolerance = 1e-6
  )

  # Test case 2
  N <- 100
  nA <- 21
  nAB <- 17
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "less", use_mid_p = TRUE, theta = theta),
    0.176809,
    tolerance = 1e-6
  )

  # Test case 3
  N <- 100
  nA <- 20
  nAB <- 10
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "less", use_mid_p = TRUE, theta = theta),
    0.0002629394,
    tolerance = 1e-6
  )
})


test_that("hwe_exact_test calculates correct p-values without mid-p correction", {
  # Test case 1 without use_mid_p
  N <- 100
  nA <- 21
  nAB <- 19
  theta <- 4
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "less", use_mid_p = FALSE, theta = theta),
    0.6903963,
    tolerance = 1e-6
  )

  # Test case 2 without use_mid_p
  N <- 100
  nA <- 21
  nAB <- 17
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "less", use_mid_p = FALSE, theta = theta),
    0.2840415,
    tolerance = 1e-6
  )

  # Test case 3 without use_mid_p
  N <- 100
  nA <- 20
  nAB <- 10
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "less", use_mid_p = FALSE, theta = theta),
    0.0005043375,
    tolerance = 1e-6
  )
})


test_that("hwe_exact_test calculates correct p-values for 'greater' alternative without mid-p correction", {
  # Test case 1 for 'greater' without use_mid_p
  N <- 100
  nA <- 21
  nAB <- 19
  theta <- 4
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "greater", use_mid_p = FALSE, theta = theta),
    0.7159585,
    tolerance = 1e-6
  )

  # Test case 2 for 'greater' without use_mid_p
  N <- 100
  nA <- 21
  nAB <- 17
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "greater", use_mid_p = FALSE, theta = theta),
    0.9304235,
    tolerance = 1e-6
  )

  # Test case 3 for 'greater' without use_mid_p
  N <- 100
  nA <- 20
  nAB <- 10
  expect_equal(
    hwe_exact_test(N, nA, nAB, alternative = "greater", use_mid_p = FALSE, theta = theta),
    0.9999785,
    tolerance = 1e-6
  )
})

