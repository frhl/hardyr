#' Calculate Log Probability of nAB Heterozygotes
#'
#' Computes the logarithm of the probability of observing a specific number of heterozygotes (nAB)
#' in a population, given the total number of individuals (N), the total number of minor alleles (nA),
#' and a theta parameter that reflects the expected heterozygosity.
#'
#' @param N The total number of diploid individuals in the population.
#' @param nA The total number of copies of the minor allele A in the population.
#' @param nAB The observed number of heterozygotes in the population.
#' @param theta The theta parameter, which reflects the relationship between the frequencies of
#' heterozygotes and homozygotes. It's used to adjust the expected probabilities under
#' the Hardy-Weinberg equilibrium or other genetic models.
#'
#' @return The log of the probability of observing exactly nAB heterozygotes given N, nA, and theta.
#'
#' @details This function first calculates the log of the numerator, which includes the term for
#' heterozygotes weighted by the theta parameter, and the factorial of the total number of individuals (N).
#' It then subtracts the log of the denominator, which includes the factorials of the number of AA homozygotes,
#' AB heterozygotes, and BB homozygotes. Finally, it subtracts the log of the normalizing constant C, calculated
#' by `calc_log_norm_c`, to obtain the final log probability.
#'
#' The function ensures that the calculated number of AA and BB homozygotes are integers, as required for
#' valid configurations in a diploid population.
#'
#' @examples
#' N <- 100
#' nA <- 21
#' nAB <- 17
#' theta <- 4
#' calc_log_p_hets(N, nA, nAB, theta)
#'
#' @export

calc_log_p_hets <- function(N, nA, nAB, theta) {
  
  # nB is the number of copies of the reference allele
  nB <- 2 * N - nA
  
  # AA and BB are homozygous genotypes
  nAA <- (nA - nAB) / 2
  nBB <- (nB - nAB) / 2
  
  # these must be integers
  stopifnot(nAA == as.integer(nAA) && nBB == as.integer(nBB))
  
  # calculate first term of equation
  log_numerator <- ((nAB/2) * log(theta) + lfactorial(N))
  log_denom <- (lfactorial(nAA) + lfactorial(nAB) + lfactorial(nBB))
  log_first_term <- log_numerator - log_denom 
  
  # get final probability
  log_c <- calc_log_norm_c(N, nA, nB, theta)
  log_probability <- log_first_term - log_c
  
  return(log_probability)
}


