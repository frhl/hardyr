#' Calculate Log of the Normalizing Constant C
#'
#' This function calculates the log of the normalizing constant C, used in the computation
#' of probabilities for genetic heterozygote configurations under a given theta parameter.
#' It applies the log-sum-exp trick to maintain numerical stability.
#'
#' @param N The total number of individuals in the population.
#' @param nA The total number of minor allele A in the population.
#' @param nB The total number of common allele B in the population. Calculated as 2*N - nA.
#' @param theta The theta parameter reflecting the relationship between heterozygotes and homozygotes frequencies.
#'
#' @return The log of the normalizing constant C, used in the computation of probabilities.
#'
#' @examples
#' N <- 100
#' nA <- 21
#' nB <- 2 * N - nA
#' theta <- 4
#' calc_log_normalizing_constant(N, nA, nB, theta)
#'
#' @export

calc_log_normalizing_constant <- function(N, nA, nB, theta) {
  
  log_C <- -Inf 
  is_odd <- as.integer(nA%%2==1)
  
  # Iterate over all possible n'AB values to compute the log of C
  for(n_prime_AB in seq(is_odd, nA, by = 2)) {
    n_prime_AA <- (nA - n_prime_AB) / 2
    n_prime_BB <- (nB - n_prime_AB) / 2
    
    # calculate numerator and denominate serperatly
    numerator <- (n_prime_AB/2) * log(theta) + lfactorial(N)
    denom <- (lfactorial(n_prime_AA) + lfactorial(n_prime_AB) + lfactorial(n_prime_BB))
    log_term <- numerator-denom
    
    # log-sum-exp trick for numerical stability
    max_log <- max(log_C, log_term)
    log_C <- max_log + log(exp(log_C - max_log) + exp(log_term - max_log))
  }
  
  return(log_C)
}
