#include <Rcpp.h>
#include "log_sum.h" // Directly or indirectly uses log_sum
#include "log_n_choose_k.h" // Directly uses log_n_choose_k
#include "probability_of_M_pairs.h"
using namespace Rcpp;

// [[Rcpp::export]]
double probability_of_M_pairs(int N, int K, int M) {
  double prob = log_n_choose_k(N, M) +
    log_n_choose_k(N - M, K - 2 * M) -
    log_n_choose_k(2 * N, K) +
    (K - 2 * M) * log(2);
  return prob;
}


