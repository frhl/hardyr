



calc_null_p <- function(N, nA, theta=4, FUN=mean){
  
  #write("This function is experimental/still under development", stderr())
  
  nAB <- nA
  nB <- 2 * N - nA
  nAA <- (nA - nAB) / 2
  nBB <- (nB - nAB) / 2
  
  # get all possible configurations
  is_odd <- as.integer(nA%%2==1)
  configurations <- seq(is_odd, nA, by=2)
  log_ps <- rep(NA, length(configurations))
  start_index <- which(configurations==nAB)
  
  # the trick here is to get the full distribution, and evaluate all
  # outcomes, and then take the cumualtive sum across.
  
  # get first probability
  log_p_start <- calc_log_p_hets(N, nA, nAB, theta)
  log_ps[start_index] <- log_p_start
  
  nABi <- nAB
  nAAi <- nAA
  nBBi <- nBB
  next_log_p <- log_p_start
  
  # go through lower
  for (i in rev(which(configurations<nAB))){
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
  
  # go through upper
  for (i in which(configurations>nAB)){
    next_log_p <-  next_log_p + log_prob_increase(nAAi, nBBi, nABi, theta)
    log_ps[i] <- next_log_p
    nABi  <- nABi + 2
    nAAi <- (nA - nABi) / 2
    nBBi <- (nB - nABi) / 2
  }
  
  # get cumulative sum of probabilities
  ps <- cumsum(exp(log_ps))
  probs <- unlist(hwe_cond_p_het(N, nA, configurations))
  sampled_ps <- sample(ps, 100, replace = TRUE,  prob = probs)
  return(FUN(sampled_ps))
}
