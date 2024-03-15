#' Sum Probability of M Pairs
#'
#' This function calculates the sum of probabilities of observing a given number of samples
#' with a biallelic variant in a gene, assessing deviation from Hardy-Weinberg equilibrium
#' at the gene level. It is particularly useful for evaluating the genetic structure
#' and variation within a population.
#'
#' @param N Numeric, the total number of individuals in the population.
#' @param K Numeric, the total number of mutated gene copies observed.
#' @param M Numeric, the number of samples that have a biallelic variant in the gene.
#' @param midp Logical, whether to use the mid-p adjustment in the probability calculation.
#'        The mid-p adjustment is used to correct for discrete probability distributions
#'        in hypothesis testing. Default is TRUE.
#'
#' @return A list containing:
#' \itemize{
#'   \item{gt}{Logical, indicating whether the function iterated over greater (TRUE) or
#'       lesser (FALSE) half of the possible M values.}
#'   \item{prob}{Numeric vector, the calculated probabilities for each M value.}
#'   \item{prob_midp}{Numeric vector or NULL, the mid-p adjusted probabilities if `midp` is TRUE,
#'       otherwise NULL.}
#' }
#'
#' @examples
#' # Calculate the sum probability of M pairs for a given N, K, and M, with mid-p adjustment
#' result <- sum_probability_of_M_pairs(N = 100, K = 20, M = 5, midp = TRUE)
#' print(result)
#'
#' @export
#'
#' @importFrom Rcpp sourceCpp
#'
#' @useDynLib yourPackageName
#' @importMethodsFrom Rcpp evalCpp

sum_probability_of_M_pairs <- function(N, K, M, midp=TRUE)
{

  # N - individuals in population
  # K - total mutated gene copies
  # M - the number of samples that have a biallelic variant in the gene
  stopifnot(is.numeric(N))
  stopifnot(is.numeric(K))
  stopifnot(is.numeric(M))
  stopifnot(is.logical(midp))
  
  if (K>(N*2)) stop("Mutated haplotypes (K) cannot exceed total number of haplotypes in population (N*2)!")
  if ((2*M)>K) stop("Bi-allelic haplotypes (2*M) cannot exceed number of mutated haplotypes (K)!" ) 
  
  probs_midp <- NULL
  
  if (M > N/2) {
    probs <- rep(0, length(seq(M+1, min(N, floor(K/2)))))
    for(i in seq(M+1, min(N, floor(K/2)))) {
      probs[i-M] <- probability_of_M_pairs(N, K, i)
    }
    
    if (midp && length(probs) > 0) {
      probs_midp <- probs
      probs_midp[1] <- probs_midp[1] + log(1/2)
    }
    
    return(list(gt=TRUE, prob=probs, prob_midp=probs_midp))

  } else {
    probs <- rep(0, M+1)
    for(i in seq(0, M)) {
      probs[i+1] <- probability_of_M_pairs(N, K, i)
      
    }
    if (midp && length(probs) > 0){
      probs_midp <- probs
      probs_midp[i+1] <- probability_of_M_pairs(N, K, i)+log(1/2)
    }
    return(list(gt=FALSE, prob=probs, prob_midp=probs_midp))
  }
}

