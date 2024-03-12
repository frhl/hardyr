# Example tests for log_prob_decrease
test_that("log_prob_decrease returns correct log probabilities", {
  expect_equal(
    log_prob_decrease(nAA = 5, nBB = 5, nAB = 10, theta = 4),
    (log(10) + log(9)) - (log(4) + log(6) + log(6)), 
    tolerance = 1e-8
  )
})

# Example tests for log_prob_increase
test_that("log_prob_increase returns correct log probabilities", {
  expect_equal(
    log_prob_increase(nAA = 5, nBB = 5, nAB = 8, theta = 4),
    log(4) + log(5) + log(5) - (log(10) + log(9)),
    tolerance = 1e-8
  )
})


