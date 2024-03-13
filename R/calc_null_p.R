

calc_null_p <- function(N, nA, theta=4, FUN=mean){
  
  nAB <- nA
  nB <- 2 * N - nA
  nAA <- (nA - nAB) / 2
  nBB <- (nB - nAB) / 2
  
  # get all possible configurations
  is_odd <- as.integer(nA%%2==1)
  configurations <- seq(is_odd, nA, by=2)
  ps <- rep(NA, length(configurations))

  # get pvalues - could do implementation much cleaner by generating
  # the full distribution and then sampling as we go along. 
  for (i in seq_along(configurations)){
    ps[i] <- hwe_exact_test(N, nA, configurations[i], theta=theta)
  }

  # sampling probabilities
  probs <- hwe_cond_p_het(N, nA, configurations)$p
  sampled_ps <- sample(ps, 100, replace = TRUE,  prob = probs)
  return(FUN(sampled_ps))
  
}

