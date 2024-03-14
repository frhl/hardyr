#' Calculate Expected P-value under the Null Hypothesis for HWE
#'
#' Computes the expected P-value under the null hypothesis of Hardy-Weinberg Equilibrium (HWE) 
#' for a given allele frequency and sample size.  
#'
#' @param N Integer, the total number of individuals in the population.
#' @param nA Integer, the total number of occurrences of one allele in the population.
#' @param theta Numeric, the scaling parameter used in calculating log probabilities of 
#'   configurations under HWE. Default is 4.
#' @param nsim Integer, the number of simulations to perform. Each simulation generates
#'   a single P-value from the null distribution. Default value should be set by the user
#'   depending on the context
#' @param use_mid_p Logical, whether to apply a mid-p correction to the computed 
#'   P-values. Default is `FALSE`.
#'
#' @details
#' The function calculates the expected P-value under the null hypothesis by first determining
#' all possible configurations of homozygous and heterozygous individuals given the total number
#' of alleles and then calculating the log probabilities for these configurations. It takes into 
#' account the non-uniform distribution of P-values under the null hypothesis due to varying 
#' allele frequencies and sample sizes. This function must be applied on EACH marker seperately. 
#'
#' @return A numeric value representing the expected P-value under the null hypothesis of HWE. 
#'
#' @examples
#' expected_p_dist <- calc_null_p(N = 100, nA = 50)
#' print(expected_p_dist) # select one to carry forward
#'
#' @export


calc_null_p <- function(N, nA, theta=4, nsim=100){
  
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
  
  # get cumulative sum of probabilities across configurations
  probs <- unlist(hwe_cond_p_het(N, nA, configurations))
  ps <- cumsum(exp(log_ps))
  sampled_ps <- sample(ps, nsim, replace = TRUE,  prob = probs)
  return(sampled_ps)
  
}


