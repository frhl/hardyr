#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double log_sum(int i) {
  if (i <= 0) {
    return 0.0;
  } else {
    double sum = 0.0;
    for (int j = 1; j <= i; ++j) {
      sum += log(j);
    }
    return sum;
  }
}

// [[Rcpp::export]]
double log_n_choose_k(int n, int k) {
  return log_sum(n) - log_sum(k) - log_sum(n - k);
}

// [[Rcpp::export]]
double probability_of_M_pairs(int N, int K, int M) {
  double prob = log_n_choose_k(N, M) +
    log_n_choose_k(N - M, K - 2 * M) -
    log_n_choose_k(2 * N, K) +
    (K - 2 * M) * log(2);
  return prob;
}

