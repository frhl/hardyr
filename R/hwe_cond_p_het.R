#' Conditional Probability of Heterozygotes under HWE
#'
#' Calculates the conditional probability of observing a specified number of
#' heterozygous individuals under Hardy-Weinberg Equilibrium (HWE) conditions
#' given the total allele counts. This function is crucial for understanding
#' genetic variation and allele distribution in populations under the assumption
#' of HWE.
#'
#' @param n Integer, the total number of individuals in the population.
#' @param nA Integer, the total number of occurrences of one allele (A) in the
#'   population. This count includes both homozygotes (AA) and heterozygotes (AB).
#' @param nAB Integer, the total number of heterozygous individuals (AB) in the
#'   population.
#'
#' @details The function first calculates the number of AA and BB homozygotes
#' based on given inputs. It then computes the probability expression's numerator
#' and denominator, utilizing factorial logarithms to accommodate large numbers.
#' The final probability is obtained by exponentiating the difference between
#' the numerator and denominator.
#'
#' The Hardy-Weinberg Equilibrium principle underpins the calculation, assuming
#' no evolutionary forces are acting on the allele frequencies (e.g., no selection,
#' mutation, migration, or genetic drift).
#'
#' @return A numeric value representing the conditional probability of observing
#' the specified number of heterozygotes (`nAB`) given the total number of
#' individuals (`n`) and the total number of one allele (`nA`), under Hardy-Weinberg
#' Equilibrium conditions.
#'
#' @examples
#' # Calculate the probability for a population of 100 individuals,
#' # with 50 occurrences of allele A and 20 heterozygotes.
#' prob <- hwe_cond_p_het(100, 50, 20)
#' print(prob)
#'
#' @note This function assumes that the population is in Hardy-Weinberg Equilibrium,
#' which may not hold in all real-world populations due to evolutionary influences.
#'
#' @references
#' For more information, see the documentation for the HardyWeinberg package on CRAN:
#' https://cran.r-project.org/web/packages/HardyWeinberg/index.html
#'
#' @export

hwe_cond_p_het <- function (n, nA, nAB) {

    nAA <- 0.5 * (nA - nAB)
    nB <- 2 * n - nA
    nBB <- 0.5 * (nB - nAB)
    numer <- lfactorial(n) + lfactorial(nA) + lfactorial(2 * 
        n - nA) + nAB * log(2)
    denom <- lfactorial(2 * n) + lfactorial(nAB) + lfactorial(0.5 * 
        (nA - nAB)) + lfactorial(0.5 * (nB - nAB))
    prob <- exp(numer - denom)
    return(prob)
}
