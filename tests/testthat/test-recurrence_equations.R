
tolerance <- 10e-10

test_that("recurrence relations produce expected values", {
  N <- 100
  nA <- 21
  nAB <- 13
  nB <- 2 * N - nA
  theta <- 2

  nAA <- (nA - nAB) / 2
  nBB <- (nB - nAB) / 2

  # get all possible configurations
  is_odd <- as.integer(nA %% 2 == 1)
  configurations <- seq(is_odd, nA, by = 2)
  log_ps <- rep(NA, length(configurations))
  start_index <- which(configurations == nAB)
  log_p_start <- calc_log_p_hets(N, nA, nAB, theta)
  log_ps[start_index] <- log_p_start

  expected_log_ps <- unlist(lapply(seq_along(configurations), function(i){
    return(calc_log_p_hets(N, nA, configurations[i], theta))
  }))

  nABi <- nAB
  nAAi <- nAA
  nBBi <- nBB
  next_log_p <- log_p_start

  # go through lower configurations
  for (i in rev(which(configurations < nAB))){
    next_log_p <- next_log_p + log_prob_decrease(nAAi, nBBi, nABi, theta)
    log_ps[i] <- next_log_p
    nABi  <- nABi - 2
    nAAi <- (nA - nABi) / 2
    nBBi <- (nB - nABi) / 2
  }

  nABi <- nAB
  nAAi <- nAA
  nBBi <- nBB
  next_log_p <- log_p_start

  # go through upper configurations
  for (i in which(configurations > nAB)){
    next_log_p <- next_log_p + log_prob_increase(nAAi, nBBi, nABi, theta)
    log_ps[i] <- next_log_p
    nABi  <- nABi + 2
    nAAi <- (nA - nABi) / 2
    nBBi <- (nB - nABi) / 2
  }

  # Check if the calculated log probabilities using recurrence are close enough to the expected ones
  differences <- abs(log_ps - expected_log_ps)
  expect_true(all(differences < tolerance), 
              info = paste("Differences are larger than the tolerance level:", tolerance))
})
