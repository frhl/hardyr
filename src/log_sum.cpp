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

