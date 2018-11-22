// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// burkardt_polynomial_dif
Rcpp::List burkardt_polynomial_dif(Rcpp::NumericVector c1, Rcpp::IntegerMatrix Powers, Rcpp::IntegerVector dif);
RcppExport SEXP _mpolynomial_burkardt_polynomial_dif(SEXP c1SEXP, SEXP PowersSEXP, SEXP difSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type Powers(PowersSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dif(difSEXP);
    rcpp_result_gen = Rcpp::wrap(burkardt_polynomial_dif(c1, Powers, dif));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_differentiate
Rcpp::List rcpp_differentiate(Rcpp::NumericVector Coeffs, Rcpp::IntegerVector Degrees, Rcpp::IntegerVector DWR);
RcppExport SEXP _mpolynomial_rcpp_differentiate(SEXP CoeffsSEXP, SEXP DegreesSEXP, SEXP DWRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Coeffs(CoeffsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Degrees(DegreesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type DWR(DWRSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_differentiate(Coeffs, Degrees, DWR));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_evaluate
Rcpp::NumericVector rcpp_evaluate(Rcpp::NumericVector Coeffs, Rcpp::IntegerVector Degrees, Rcpp::NumericMatrix Values);
RcppExport SEXP _mpolynomial_rcpp_evaluate(SEXP CoeffsSEXP, SEXP DegreesSEXP, SEXP ValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Coeffs(CoeffsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Degrees(DegreesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Values(ValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_evaluate(Coeffs, Degrees, Values));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_rank_grlex
int rcpp_rank_grlex(Rcpp::NumericVector Powers);
RcppExport SEXP _mpolynomial_rcpp_rank_grlex(SEXP PowersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Powers(PowersSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_rank_grlex(Powers));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_unrank_grlex
Rcpp::IntegerVector rcpp_unrank_grlex(int m, int rank);
RcppExport SEXP _mpolynomial_rcpp_unrank_grlex(SEXP mSEXP, SEXP rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type rank(rankSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_unrank_grlex(m, rank));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mpolynomial_burkardt_polynomial_dif", (DL_FUNC) &_mpolynomial_burkardt_polynomial_dif, 3},
    {"_mpolynomial_rcpp_differentiate", (DL_FUNC) &_mpolynomial_rcpp_differentiate, 3},
    {"_mpolynomial_rcpp_evaluate", (DL_FUNC) &_mpolynomial_rcpp_evaluate, 3},
    {"_mpolynomial_rcpp_rank_grlex", (DL_FUNC) &_mpolynomial_rcpp_rank_grlex, 1},
    {"_mpolynomial_rcpp_unrank_grlex", (DL_FUNC) &_mpolynomial_rcpp_unrank_grlex, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_mpolynomial(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}