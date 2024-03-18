#' Perform an Exact Test under Hardy-Weinberg Equilibrium
#'
#' Conducts an exact test for Hardy-Weinberg equilibrium for a given set of genotypic frequencies,
#' specified by the number of individuals (N), the number of minor alleles (nA), and the observed
#' number of heterozygotes (nAB), considering a specified inbreeding coefficient (theta).
#' The function allows for testing whether the observed heterozygosity is significantly less than,
#' greater than, or different from the expected heterozygosity under Hardy-Weinberg equilibrium.
#'
#' @param N Numeric, the total number of individuals in the population.
#' @param K Numeric, the total number of mutated gene copies observed (equivalent to minor allele count)
#' @param M Numeric, the number of samples that have are homozygous carriers for the minor allele
#' @param theta Numeric, the inbreeding coefficient parameter; defaults to 4. This parameter reflects 
#' the degree of deviation from random mating.
#' @param alternative Character, specifies the alternative hypothesis and must be one of "less", 
#' "greater", or "two.sided". "Less" tests for a deficit of heterozygotes (indicating excess homozygosity),
#' and "greater" for an excess of heterozygotes.
#' @param use_mid_p Logical, indicating whether to use the mid-p correction for the p-value calculation in one-sided tests; defaults to TRUE.
#' @param debug Logical, if TRUE, the function will enter debug mode for troubleshooting; defaults to FALSE.
#' 
#' @return Numeric, the p-value for the test.
#'
#' @examples
#' N <- 100
#' K <- 21 # minor allele count
#' M <- 3 # 3 homozygous carriers
#' p_value <- hwe_exact_test(N, K, M, theta = 4, alternative = "less")
#'
#' @export


hwe_exact_test <- function(N, K, M, theta=4, alternative="less", use_mid_p=TRUE, debug=FALSE){
  
  if (debug) browser()
  
  # some tests
  allele_flip <- FALSE
  stopifnot(alternative %in% c("less", "greater"))
  stopifnot(is.numeric(N) && N > 0)
  stopifnot(is.numeric(K) && K >= 0)
  stopifnot(is.numeric(M) && M >= 0)
  stopifnot(is.numeric(theta))
  stopifnot(is.logical(use_mid_p))
  stopifnot(is.logical(debug))
  
  # ensure that counts line up
  if (K>(N*2)) stop("Mutated haplotypes (K) cannot exceed total number of haplotypes in population (N*2)!")
  if ((2*M)>K) stop("Bi-allelic haplotypes (2*M) cannot exceed number of mutated haplotypes (K)!" ) 
  
  # check if A is minor alleles otherwise flip
  actual_mac <- get_mac(N, K, M)
  if (actual_mac != K){
    allele_flip <- TRUE
    allele_counts <- get_allele_counts(N, K, M)
    K <- get_mac(N, K, M)
    M <- allele_counts[3]
    warning(paste0("param 'K' is not the minor allele! Using K=",K," and M=",M, " instead!"))
  }
  
  # map to nA, nAB, nBB
  nAB <- K - M*2
  nA <- K
  
  # always estimate these guys
  nB <- 2 * N - nA
  nAA <- (nA - nAB) / 2
  nBB <- (nB - nAB) / 2
  
  # check allele counts
  stopifnot(nAA + nAB + nBB == N)
  stopifnot(nAA>=0, nBB>=0, nB>=0)
  
  # get all possible configurations
  is_odd <- as.integer(nA%%2==1)
  configurations <- seq(is_odd, nA, by=2)
  log_ps <- rep(NA, length(configurations))
  start_index <- which(configurations==nAB)
  
  # get first probability
  log_p_start <- calc_log_p_hets(N, nA, nAB, theta)
  log_ps[start_index] <- log_p_start
  
  # use recurrent relation to calculate lower
  if (alternative %in% c("less", "two.sided") ){
    
    nABi <- nAB
    nAAi <- nAA
    nBBi <- nBB
    next_log_p <- log_p_start
    
    # go through lower
    for (i in rev(which(configurations<nAB))){
      next_log_p <-  next_log_p + log_prob_decrease(nAAi, nBBi, nABi, theta)
      log_ps[i] <- next_log_p
      nABi  <- nABi -2
      nAAi <- (nA - nABi) / 2
      nBBi <- (nB - nABi) / 2
    }
  }
  
  # use recurrent relation to calculate upper
  if (alternative %in% c("greater", "two.sided")){
    
    nABi <- nAB
    nAAi <- nAA
    nBBi <- nBB
    next_log_p <- log_p_start
    
    # go through upper
    for (i in which(configurations>nAB)){
      next_log_p <-  next_log_p + log_prob_increase(nAAi, nBBi, nABi, theta)
      log_ps[i] <- next_log_p
      nABi  <- nABi + 2
      nAAi <- (nA - nABi) / 2
      nBBi <- (nB - nABi) / 2
    }
  }
  
  # get mid P value for one-sided test
  if (use_mid_p) log_ps[start_index] <- log(1/2) + log_ps[start_index]
  
  # deal with sides
  if (alternative %in% "less"){
    p <- sum(exp(log_ps[which(configurations<=nAB)]))
  } else if (alternative %in% "greater"){
    p <- sum(exp(log_ps[which(configurations>=nAB)]))
  }
  
  return(p)
}
