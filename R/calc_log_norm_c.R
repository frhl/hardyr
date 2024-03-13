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
#' calc_log_norm_c(N, nA, nB, theta)
#'
#' @export

calc_log_norm_c <- function(N, nA, nB, theta) {
  
  log_C <- -Inf 
  is_odd <- as.integer(nA%%2==1)
  log_terms <- c()
  # Iterate over all possible n'AB values to compute the log of C
  for(n_prime_AB in seq(is_odd, nA, by = 2)) {
    n_prime_AA <- (nA - n_prime_AB) / 2
    n_prime_BB <- (nB - n_prime_AB) / 2
    
    # calculate numerator and denominate serperatly
    numerator <- (n_prime_AB/2) * log(theta) + lfactorial(N)
    denom <- (lfactorial(n_prime_AA) + lfactorial(n_prime_AB) + lfactorial(n_prime_BB))
    value <- numerator-denom
    log_terms <- c(log_terms, value)
  }
  log_C <- log_sum_exp(log_terms)
  return(log_C)
}



