#' Calculate Expected P-values from Uniform
#'
#' This function computes the expected uniform p-values based on the observed p-values. 
#' It can optionally remove NA values before computation. The function is designed 
#' to transform a set of observed p-values such that they follow a uniform distribution, 
#' which is useful for certain statistical analyses where uniformity of p-values is expected under the null hypothesis.
#'
#' @param observed_p A numeric vector of observed p-values. Must be non-negative.
#' @param na.rm Logical; if TRUE, any NA values in `observed_p` will be removed before computation.
#'               Default is FALSE, which will issue a warning if NAs are present but not remove them.
#'
#' @details The function first checks that `observed_p` is a numeric vector and contains 
#'          non-negative values. It then handles NA values according to the `na.rm` parameter. 
#'          The expected uniform p-values are calculated by ranking the observed p-values and 
#'          then scaling these ranks linearly between 1/n and 1, where n is the number of 
#'          observed p-values after NA handling. This transformation assumes that the input p-values 
#'          are drawn from a continuous uniform distribution.
#'
#' @return Returns a numeric vector of the same length as `observed_p` (after NA removal, if applicable), 
#'         containing the expected uniform p-values corresponding to the ranks of the observed p-values.
#'
#' @examples
#' # Example without NAs and na.rm = FALSE
#' observed_p <- c(0.05, 0.2, 0.01, 0.5, 0.95)
#' calc_expected_p_uniform(observed_p)
#'
#' # Example with NAs and na.rm = TRUE
#' observed_p_with_na <- c(0.05, NA, 0.01, 0.5, 0.95)
#' calc_expected_p_uniform(observed_p_with_na, na.rm = TRUE)
#'
#' @export
calc_expected_p_uniform <- function(observed_p, na.rm = FALSE){
  stopifnot(is.numeric(observed_p))
  stopifnot(any(observed_p >= 0))
  sum_na <- sum(is.na(observed_p))
  if (sum_na > 0){
    if (na.rm){
      observed_p <- observed_p[!is.na(observed_p)]
    } else {
      warning("NAs detected in P-values. Set na.rm = TRUE to remove them.")
    }
  }
  n <- length(observed_p)
  observed_rank <- rank(observed_p)
  uniform <- (1:n)/(n+1)
  uniform <- uniform[observed_rank]
  return(uniform)
}

