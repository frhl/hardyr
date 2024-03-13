#' Calculate Theta using Inbreeding Coeffecient
#'
#' Computes the theta parameter used in Hardy-Weinberg equilibrium tests,
#' taking into account the inbreeding coefficient.
#'
#' @param pA Numeric, the frequency of the minor allele in the population.
#' @param f Numeric, the inbreeding coefficient.
#'
#' @return Numeric, the calculated theta value.
#'
#' @examples
#' pA <- 0.2  # Minor allele frequency
#' f <- 0.01 # Inbreeding coefficient
#' theta <- calc_theta_from_f(pA, f)
#' print(theta)
#'
#' @export
calc_theta_from_f <- function(pA, f) {
  stopifnot(is.numeric(pA), is.numeric(f), pA >= 0, pA <= 1)
  num <- 4 * pA * (1 - pA) * (1 - f)^2
  den <- (pA + f * (1 - pA)) * ((1 - pA) + f * pA)
  theta <- num/den
  return(theta)
}


