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
// burkardt_polynomial_value
Rcpp::NumericVector burkardt_polynomial_value(Rcpp::NumericVector c, Rcpp::IntegerMatrix Powers, Rcpp::NumericMatrix x);
RcppExport SEXP _mpolynomial_burkardt_polynomial_value(SEXP cSEXP, SEXP PowersSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type Powers(PowersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(burkardt_polynomial_value(c, Powers, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mpolynomial_burkardt_polynomial_dif", (DL_FUNC) &_mpolynomial_burkardt_polynomial_dif, 3},
    {"_mpolynomial_burkardt_polynomial_value", (DL_FUNC) &_mpolynomial_burkardt_polynomial_value, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_mpolynomial(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
