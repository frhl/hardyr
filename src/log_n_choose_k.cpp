#include <Rcpp.h>
#include "log_sum.h"
#include "log_n_choose_k.h"

using namespace Rcpp;

// [[Rcpp::export]]
double log_n_choose_k(int n, int k) {
  return log_sum(n) - log_sum(k) - log_sum(n - k);
}
