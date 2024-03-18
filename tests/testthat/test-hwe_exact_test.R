test_that("hwe_exact_test calculates correct p-values", {
  theta <- 4

  # Converted Test case 1
  N <- 100
  K <- 21  # Equivalent to nA
  M <- (21 - 19) / 2  # Calculated from nAB = K - M*2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "less", use_mid_p = TRUE),
    0.4872189,
    tolerance = 1e-6
  )

  # Converted Test case 2
  N <- 100
  K <- 21
  M <- (21 - 17) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "less", use_mid_p = TRUE),
    0.176809,
    tolerance = 1e-6
  )

  # Converted Test case 3
  N <- 100
  K <- 20
  M <- (20 - 10) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "less", use_mid_p = TRUE),
    0.0002629394,
    tolerance = 1e-6
  )
})

test_that("hwe_exact_test calculates correct p-values without mid-p correction", {
  theta <- 4

  # Converted Test case 1 without use_mid_p
  N <- 100
  K <- 21
  M <- (21 - 19) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "less", use_mid_p = FALSE),
    0.6903963,
    tolerance = 1e-6
  )

  # Converted Test case 2 without use_mid_p
  N <- 100
  K <- 21
  M <- (21 - 17) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "less", use_mid_p = FALSE),
    0.2840415,
    tolerance = 1e-6
  )

  # Converted Test case 3 without use_mid_p
  N <- 100
  K <- 20
  M <- (20 - 10) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "less", use_mid_p = FALSE),
    0.0005043375,
    tolerance = 1e-6
  )
})

test_that("hwe_exact_test calculates correct p-values for 'greater' alternative without mid-p correction", {
  theta <- 4

  # Converted Test case 1 for 'greater' without use_mid_p
  N <- 100
  K <- 21
  M <- (21 - 19) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "greater", use_mid_p = FALSE),
    0.7159585,
    tolerance = 1e-6
  )

  # Converted Test case 2 for 'greater' without use_mid_p
  N <- 100
  K <- 21
  M <- (21 - 17) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "greater", use_mid_p = FALSE),
    0.9304235,
    tolerance = 1e-6
  )

  # Converted Test case 3 for 'greater' without use_mid_p
  N <- 100
  K <- 20
  M <- (20 - 10) / 2
  expect_equal(
    hwe_exact_test(N, K, M, theta = theta, alternative = "greater", use_mid_p = FALSE),
    0.9999785,
    tolerance = 1e-6
  )
})




test_that("hwe_exact_test output matches log_sum_probability_of_M_pairs for odd values of K", {
  N = 100; K=21; M=0;
  expect_equal(
    hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE),
    sum(exp(log_sum_probability_of_M_pairs(N = N, K = K, M = M, midp = TRUE)$prob_midp)),
    tolerance = 1e-9
  )
  
  N = 100; K=21; M=1;
  expect_equal(
    hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE),
    sum(exp(log_sum_probability_of_M_pairs(N = N, K = K, M = M, midp = TRUE)$prob_midp)),
    tolerance = 1e-9
  )
  
  N = 100; K=21; M=2;
  expect_equal(
    hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE),
    sum(exp(log_sum_probability_of_M_pairs(N = N, K = K, M = M, midp = TRUE)$prob_midp)),
    tolerance = 1e-9
  )
  
  N = 100; K=21; M=3;
  expect_equal(
    hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE),
    sum(exp(log_sum_probability_of_M_pairs(N = N, K = K, M = M, midp = TRUE)$prob_midp)),
    tolerance = 1e-9
  )
  
  N = 100; K=21; M=10;
  expect_equal(
    hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE),
    sum(exp(log_sum_probability_of_M_pairs(N = N, K = K, M = M, midp = TRUE)$prob_midp)),
    tolerance = 1e-9
  )
  
  N = 100; K=21; M=11;
  expect_error(
    hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE)
  )

})

test_that("hwe_exact_test output matches log_sum_probability_of_M_pairs for even values of K", {
  N = 100; K = 22; # Using an even value for K
  M_values = c(0, 1, 2, 3, 10, 11) # The same range of M values
  
  for(M in M_values) {
    if(2 * M <= K) {
      expect_equal(
        hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE),
        sum(exp(log_sum_probability_of_M_pairs(N = N, K = K, M = M, midp = TRUE)$prob_midp)),
        tolerance = 1e-9
      )
    } else {
      expect_error(
        hwe_exact_test(N = N, K = K, M = M, theta = 4, alternative = "greater", use_mid_p = TRUE)
      )
    }
  }
})

