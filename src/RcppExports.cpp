// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// IT
IntegerMatrix IT(NumericVector rowsums, NumericVector colsums, int nrow, int ncol, int N);
RcppExport SEXP _spasm_IT(SEXP rowsumsSEXP, SEXP colsumsSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rowsums(rowsumsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colsums(colsumsSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(IT(rowsums, colsums, nrow, ncol, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spasm_IT", (DL_FUNC) &_spasm_IT, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_spasm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
