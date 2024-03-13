#' Calculate Normalizing Constant
#'
#' Computes the normalizing constant \( C \) given the population size \( N \), the number of alleles \( nA \), \( nB \), and inbreeding coefficient \( theta \).
#'
#' @param N The total number of individuals in the population.
#' @param nA The number of alleles A in the population.
#' @param nB The number of alleles B in the population.
#' @param theta The inbreeding coefficient.
#'
#' @return The normalizing constant \( C \).
#' @export
#' @examples
#' calc_norm_c(100, 40, 160, 0.01)


calc_norm_c <- function(N, nA, nB, theta) {
  
  terms <- c()
  is_odd <- as.integer(nA%%2==1)
  for(n_prime_AB in seq(is_odd, nA, by = 2)) {
    n_prime_AA <- (nA - n_prime_AB) / 2
    n_prime_BB <- (nB - n_prime_AB) / 2
    numerator <- (theta^(n_prime_AB/2))*factorial(N)
    denom <- (factorial(n_prime_AA) * factorial(n_prime_AB) * factorial(n_prime_BB))
    value <- numerator/denom
    terms <- c(terms, value)
  }
  constant <- sum(terms)
  return(constant)
}
