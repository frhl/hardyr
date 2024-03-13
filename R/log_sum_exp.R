#' Log-Sum-Exp Trick
#'
#' Performs stable computation of the logarithm of the sum of exponentials of input log values.
#' This function is particularly useful when dealing with log probabilities to prevent
#' numerical underflow or overflow. It also handles cases where input values may include `-Inf`.
#'
#' @param log_values Numeric vector, log-transformed values to be summed.
#'
#' @return Numeric, the log of the sum of exponentials of the input values.
#'
#' @examples
#' log_probs <- c(-Inf, log(0.01), log(0.99))
#' log_sum <- log_sum_exp(log_probs)
#' print(log_sum)
#'
#' @export

log_sum_exp <- function(log_values) {
  # Filter out non-finite values to prevent issues with `-Inf`
  finite_log_values <- log_values[is.finite(log_values)]
  if (length(finite_log_values) == 0) {
    return(-Inf)
  }
  
  # Find the maximum log value to scale the computation and avoid overflow
  max_log <- max(finite_log_values)
  # Scale the log values by subtracting the maximum log value and exponentiate
  scaled_log_values <- exp(finite_log_values - max_log)
  
  # Sum the scaled values and take the log, then re-scale by adding max_log
  return(max_log + log(sum(scaled_log_values)))
}

