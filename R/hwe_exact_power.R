#' Calculate the power of an exact test for Hardy-Weinberg Equilibrium
#'
#' Computes the statistical power of the exact test for Hardy-Weinberg Equilibrium
#' given a specific number of individuals, the number of minor alleles, and a theta parameter.
#' The power is defined as the probability that the test correctly rejects the null hypothesis
#' (HWE) when the alternative hypothesis is true. For a given θ, the power is calculated by summing 
#' the probabilities of all genotype counts that lead to the rejection of the null hypothesis.
#' While the rejection region is determined under the null hypothesis, the power calculation
#' is performed under the alternative hypothesis conditions.
#'
#' @param N Integer, total number of diploid individuals in the population.
#' @param K Integer, total number of copies of the minor allele in the population (equivalent to minor allele count)
#' @param theta Numeric, optional, deviation parameter reflecting inbreeding coefficient; if not provided, `f` must be specified.
#' @param f Numeric, optional, inbreeding coefficient; if not provided, `theta` must be specified. Used to calculate `theta` if `theta` is not directly provided.
#' @param sig.level Numeric, significance level for identifying the rejection region under the null hypothesis; defaults to 0.05.
#' @param use_mid_p Logical, indicates whether to use the mid-p correction for the p-value calculation; defaults to FALSE. The mid-p correction is applied only to HWE p-values.
#' @param alternative Character, specifies the direction of the alternative hypothesis; currently only "less" and "greater" is implemented.
#'
#'
#' @return Numeric, the power of the Hardy-Weinberg Equilibrium exact test.
#'
#' @details
#' The power calculation involves:
#' 1. Calculating log probabilities of genotype configurations under both the null hypothesis (theta = 4)
#'    and the specified alternative hypothesis (current theta).
#' 2. Determining the rejection region based on probabilities under the null hypothesis.
#' 3. Summing probabilities under the alternative hypothesis for configurations within the rejection region.
#'
#' This methodology allows assessing how effectively the HWE test can detect deviations
#' from equilibrium for different theta values, representing various inbreeding coefficients
#' or other forms of non-random mating.
#'
#' @examples
#' # Calculate power with specified theta
#' power <- hwe_exact_power(N = 100, K = 14, theta = 2, sig.level = 0.05, alternative="greater")
#'
#' # Calculate power with specified inbreeding coefficient f
#' power_f <- hwe_exact_power(N = 100, K = 14, f = 0.01, sig.level = 0.05)
#'
#' @export

hwe_exact_power <- function(N, K, theta=NULL, f=NULL, sig.level=0.05, alternative="less", use_mid_p=FALSE){

  stopifnot(alternative %in% c("less", "greater"))
  if (!is.null(theta) & !is.null(f)) stop("Only one of param 'theta' or 'f' can be set!")
  if (K>N) stop("'K' is not the minor allele count!")
  
  # map to nA, nAB, nBB
  nA <- K
  nAB <- nA
  nB <- 2 * N - nA
  nAA <- (nA - nAB) / 2
  nBB <- (nB - nAB) / 2
  
  # calculate theta if not specified
  if ((is.null(theta) & (!is.null(f)))){
    pA <- nA / (2 * N)
    theta <- calc_theta_from_f(pA, f)
  }

  # get all possible configurations
  is_odd <- as.integer(nA%%2==1)
  configurations <- seq(is_odd, nA, by=2)
  log_ps <- rep(NA, length(configurations))
  log_ps_hwe <- rep(NA, length(configurations))
  start_index <- which(configurations==nAB)
  
  # get first probability
  log_p_start <- calc_log_p_hets(N, nA, nAB, theta)
  log_ps[start_index] <- log_p_start
  
  # Also calculate the first probability when theta=4
  log_p_start_hwe <- calc_log_p_hets(N, nA, nAB, theta=4)
  log_ps_hwe[start_index] <- log_p_start_hwe
  
  nABi <- nAB
  nAAi <- nAA
  nBBi <- nBB
  next_log_p <- log_p_start
  next_log_p_hwe <- log_p_start_hwe
  
  # go through lower
  for (i in rev(which(configurations<nAB))){
    next_log_p <- next_log_p + log_prob_decrease(nAAi, nBBi, nABi, theta)
    next_log_p_hwe <- next_log_p_hwe + log_prob_decrease(nAAi, nBBi, nABi, theta=4)
    log_ps[i] <- next_log_p
    log_ps_hwe[i] <- next_log_p_hwe
    nABi  <- nABi - 2
    nAAi <- (nA - nABi) / 2
    nBBi <- (nB - nABi) / 2
  }

  nABi <- nAB
  nAAi <- nAA
  nBBi <- nBB
  next_log_p <- log_p_start
  next_log_p_hwe <- log_p_start_hwe
  
  # go through upper
  for (i in which(configurations>nAB)){
    next_log_p <-  next_log_p + log_prob_increase(nAAi, nBBi, nABi, theta)
    next_log_p_hwe <-  next_log_p_hwe + log_prob_increase(nAAi, nBBi, nABi, theta=4)
    log_ps[i] <- next_log_p
    log_ps_hwe[i] <- next_log_p_hwe
    nABi  <- nABi + 2
    nAAi <- (nA - nABi) / 2
    nBBi <- (nB - nABi) / 2
  }
  
  # mid p correction only for hwe p values
  if (use_mid_p) {
    log_ps_hwe[start_index] <- log(1/2) + log_ps_hwe[start_index]
  }
  
  # Determine rejection criteria based on HWE (theta = 4)
  # Calculate power by summing probabilities under the 
  # current theta for rejected configurations
  if (alternative == "greater"){
     ps_hwe_cumsum <- rev(cumsum(rev(exp(log_ps_hwe))))
  } else if (alternative == "less"){
     ps_hwe_cumsum <- cumsum(exp(log_ps_hwe))
  }
  
  rejection_criteria <- which(ps_hwe_cumsum < sig.level)
  power <- sum(exp(log_ps[rejection_criteria]))
  
  return(power)
}

