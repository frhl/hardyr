


# Test for calc_norm_c
test_that("calc_norm_c returns accurate non-log values", {
  N <- 100
  nA <- 8
  nB <- 32
  thetas <- c(0.01, 0.1, 2, 0.5, 1, 4, 8, 10)

  for (theta in thetas) {
    log_result <- exp(calc_log_norm_c(N, nA, nB, theta))
    non_log_result <- calc_norm_c(N, nA, nB, theta)
    expect_equal(log_result, non_log_result, tolerance = .Machine$double.eps^0.5)
  }

  # Larger values of N
  N <- 200
  for (theta in thetas) {
    log_result <- exp(calc_log_norm_c(N, nA, nB, theta))
    non_log_result <- calc_norm_c(N, nA, nB, theta)
    expect_equal(log_result, non_log_result, tolerance = 0.01)  # Increased tolerance for larger numbers
  }
})


# Test for calc_norm_c
test_that("calc_norm_c returns accurate non-log values", {
  N <- 100
  nA <- 14
  nB <- 2 * N - nA
  thetas <- c(0.01, 0.1, 2, 4, 6, 8, 10)
  for (theta in thetas) {
    log_result <- exp(calc_log_norm_c(N, nA, nB, theta))
    non_log_result <- calc_norm_c(N, nA, nB, theta)
    expect_equal(log_result, non_log_result, tolerance = .Machine$double.eps^0.5)
  }
})

# Test for calc_norm_c at extreme cases
test_that("calc_norm_c returns accurate non-log values when N=nA=100", {
  N <- 100
  nA <- 100
  nB <- 2 * N - nA
  thetas <- c(0.01, 0.1, 2, 4, 6, 8, 10)
  for (theta in thetas) {
    log_result <- exp(calc_log_norm_c(N, nA, nB, theta))
    non_log_result <- calc_norm_c(N, nA, nB, theta)
    expect_equal(log_result, non_log_result, tolerance = .Machine$double.eps^0.5)
  }
})




