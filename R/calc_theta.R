#' Calculate Theta
#'
#' This function calculates the theta value based on allele counts. Theta is calculated using
#' the formula \eqn{\theta = \frac{p_{AB}^2}{p_{AA} \cdot p_{BB}}}, where \eqn{p_{AA}}, \eqn{p_{AB}}, and
#' \eqn{p_{BB}} are the frequencies of the homozygous for the major allele, heterozygous, and homozygous for the
#' minor allele genotypes, respectively. These frequencies are derived from allele counts obtained
#' for a given population.
#'
#' @param N The total number of samples in the population. This should be a positive integer.
#' @param K The count of the minor allele in the population. This should be a positive integer less than or equal to 2*N.
#' @param M The number of homozygous carriers of the minor allele in the population. This should be a positive integer less than or equal to N.
#'
#' @return Returns the theta value as a numeric value. This value is a measure of the balance between
#'         the frequencies of the homozygous major, heterozygous, and homozygous minor genotypes.
#'
#' @examples
#' N <- 10000 # samples in population
#' K <- 1000 # minor allele count
#' M <- 10 # homozygous carriers
#' theta <- calc_theta(N, K, M)
#'
#' @export

calc_theta <- function(N, K, M){
  allele_counts <- get_allele_counts(N, K, M)
  n <- sum(allele_counts)
  pAA <- allele_counts[1]/n
  pAB <- allele_counts[2]/n
  pBB <- allele_counts[3]/n
  theta <- (pAB^2)/(pAA*pBB)
  return(theta)
}


