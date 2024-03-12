#' Calculate Log Probability Decrease for Heterozygotes
#'
#' Computes the decrease in log probability when the number of heterozygotes decreases by two,
#' given the counts of homozygous minor (nAA), homozygous major (nBB) genotypes, current heterozygotes (nAB),
#' and the theta parameter.
#'
#' @param nAA Integer, the number of homozygous minor genotypes.
#' @param nBB Integer, the number of homozygous major genotypes.
#' @param nAB Integer, the current number of heterozygotes.
#' @param theta Numeric, the theta parameter reflecting the relationship between heterozygotes and homozygotes frequencies.
#'
#' @return Numeric, the change in log probability resulting from a decrease of two heterozygotes.

log_prob_decrease <- function(nAA, nBB, nAB, theta){
  stopifnot(nAB >= 0)
  return((log(nAB)+log(nAB-1))-(log(theta)+log(nAA+1)+log(nBB+1)))
}

#' Calculate Log Probability Increase for Heterozygotes
#'
#' Computes the increase in log probability when the number of heterozygotes increases by two,
#' given the counts of homozygous minor (nAA), homozygous major (nBB) genotypes, current heterozygotes (nAB),
#' and the theta parameter.
#'
#' @param nAA Integer, the number of homozygous minor genotypes.
#' @param nBB Integer, the number of homozygous major genotypes.
#' @param nAB Integer, the current number of heterozygotes.
#' @param theta Numeric, the theta parameter reflecting the relationship between heterozygotes and homozygotes frequencies.
#'
#' @return Numeric, the change in log probability resulting from an increase of two heterozygotes.

log_prob_increase <- function(nAA, nBB, nAB, theta){
  stopifnot(nAB >= 0)
  return(log(theta)+log(nAA)+log(nBB)-(log(nAB+2)+log(nAB+1)))
}

