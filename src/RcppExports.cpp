// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_n_choose_k
double log_n_choose_k(int n, int k);
RcppExport SEXP _hardyr_log_n_choose_k(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(log_n_choose_k(n, k));
    return rcpp_result_gen;
END_RCPP
}
// log_sum
double log_sum(int i);
RcppExport SEXP _hardyr_log_sum(SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum(i));
    return rcpp_result_gen;
END_RCPP
}
// probability_of_M_pairs
double probability_of_M_pairs(int N, int K, int M);
RcppExport SEXP _hardyr_probability_of_M_pairs(SEXP NSEXP, SEXP KSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(probability_of_M_pairs(N, K, M));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hardyr_log_n_choose_k", (DL_FUNC) &_hardyr_log_n_choose_k, 2},
    {"_hardyr_log_sum", (DL_FUNC) &_hardyr_log_sum, 1},
    {"_hardyr_probability_of_M_pairs", (DL_FUNC) &_hardyr_probability_of_M_pairs, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_hardyr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
