

# this function is taken from https://cran.r-project.org/web/packages/HardyWeinberg/index.html 

hwe_cond_p_het <- function (n, nA, nAB) {

    nAA <- 0.5 * (nA - nAB)
    nB <- 2 * n - nA
    nBB <- 0.5 * (nB - nAB)
    numer <- lfactorial(n) + lfactorial(nA) + lfactorial(2 * 
        n - nA) + nAB * log(2)
    denom <- lfactorial(2 * n) + lfactorial(nAB) + lfactorial(0.5 * 
        (nA - nAB)) + lfactorial(0.5 * (nB - nAB))
    prob <- exp(numer - denom)
    return(list(p = prob))
}
