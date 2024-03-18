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


#' Get Allele Counts
#'
#' This function calculates the allele counts based on the total number of individuals in the population (N), 
#' the total number of mutated gene copies observed (K), and the number of samples that are homozygous 
#' carriers for the minor allele (M). It returns a vector with the counts of homozygous alternate, 
#' heterozygous, and homozygous reference allele carriers, respectively.
#'
#' @param N Numeric, the total number of individuals in the population.
#' @param K Numeric, the total number of mutated gene copies observed (equivalent to minor allele count).
#' @param M Numeric, the number of samples that are homozygous carriers for the minor allele.
#'
#' @export
get_allele_counts<- function(N, K, M){
  stopifnot(is.numeric(N) && N > 0)
  stopifnot(is.numeric(K) && K >= 0)
  stopifnot(is.numeric(M) && M >= 0)
  if (K>(N*2)) stop("Mutated haplotypes (K) cannot exceed total number of haplotypes in population (N*2)!")
  if ((2*M)>K) stop("Bi-allelic haplotypes (2*M) cannot exceed number of mutated haplotypes (K)!" ) 
  if ((K-(2*M))>N) stop(paste0("There are more heterozygous carriers (K-2*M=",K-(2*M),") than samples (N=",N,")!"))
  return(c(M, (K-2*M), N-(K-2*M)-M))
}


#' Get Minimum Allele Count (MAC)
#'
#' This function calculates the minimum allele count (MAC) for a given set of genotype counts. It uses
#' the total number of individuals (N), the total number of mutated gene copies observed (K), and the number
#' of samples that are homozygous carriers for the minor allele (M) to determine the MAC.
#'
#' @param N Numeric, the total number of individuals in the population.
#' @param K Numeric, the total number of mutated gene copies observed (equivalent to minor allele count).
#' @param M Numeric, the number of samples that are homozygous carriers for the minor allele.
#'
#' @export
get_mac <- function(N, K, M){
  stopifnot(is.numeric(N) && N > 0)
  stopifnot(is.numeric(K) && K >= 0)
  stopifnot(is.numeric(M) && M >= 0)
  if (K>(N*2)) stop("Mutated haplotypes (K) cannot exceed total number of haplotypes in population (N*2)!")
  if ((2*M)>K) stop("Bi-allelic haplotypes (2*M) cannot exceed number of mutated haplotypes (K)!" ) 
  if ((K-(2*M))>N) stop(paste0("There are more heterozygous carriers (K-2*M=",K-(2*M),") than samples (N=",N,")!"))
  counts <- get_allele_counts(N, K, M)
  hom_alt <- counts[1]*2
  hets <- counts[2]
  hom_ref <- counts[3]*2
  return(min(hom_alt+hets, hom_ref+hets))
}


