#' Calculate Inbreeding Coefficient from Theta
#'
#' Computes the inbreeding coefficient `f` based on the theta parameter used in
#' Hardy-Weinberg equilibrium tests and allele frequency `p`.
#'
#' @param p Numeric, the frequency of the minor allele in the population.
#' @param f Numeric, the inbreeding coefficient.
#'
#' @return Numeric, the inbreeding coefficient `f` recalculated from `theta` and `pA`.
#'
#' @examples
#' pA <- 0.2  # Minor allele frequency
#' f <- 0.01 
#' calc_f_from_theta(pA, f)
#'
#' @export
calc_f_from_theta <- function(pA, theta) {
  stopifnot(is.numeric(pA), is.numeric(theta), pA >= 0, pA <= 1)
  part <- pA * (1 - pA) * (4 - theta)
  z1 <- part
  z2 <- -2 * part - theta
  z3 <- part
  z <- c(z1, z2, z3)
  out <- polyroot(z)
  f <- Re(out)
  minoraf <- min(pA, 1 - pA)
  fmin <- -minoraf/(1 - minoraf)
  f <- f[f > fmin]
  f <- f[f < 1]
  return(f)
}
