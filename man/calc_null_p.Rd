% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_null_p.R
\name{calc_null_p}
\alias{calc_null_p}
\title{Calculate Expected P-value under the Null Hypothesis for HWE}
\usage{
calc_null_p(N, nA, theta = 4, FUN = mean)
}
\arguments{
\item{N}{Integer, the total number of individuals in the population.}

\item{nA}{Integer, the total number of occurrences of one allele in the population.}

\item{theta}{Numeric, the scaling parameter used in calculating log probabilities of 
configurations under HWE. Default is 4.}

\item{nsim}{Integer, the number of simulations to perform. Each simulation generates
a single P-value from the null distribution. Default value should be set by the user
depending on the context}

\item{use_mid_p}{Logical, whether to apply a mid-p correction to the computed 
P-values. Default is `FALSE`.}
}
\value{
A numeric value representing the expected P-value under the null hypothesis of HWE.
}
\description{
Computes the expected P-value under the null hypothesis of Hardy-Weinberg Equilibrium (HWE) 
for a given allele frequency and sample size.
}
\details{
The function calculates the expected P-value under the null hypothesis by first determining
all possible configurations of homozygous and heterozygous individuals given the total number
of alleles and then calculating the log probabilities for these configurations. It takes into 
account the non-uniform distribution of P-values under the null hypothesis due to varying 
allele frequencies and sample sizes. This function must be applied on EACH marker seperately.
}
\examples{
expected_p_dist <- calc_null_p(N = 100, nA = 50)
print(expected_p_dist) # select one to carry forward

}
