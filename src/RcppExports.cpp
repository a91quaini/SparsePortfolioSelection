// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// add_ints
int add_ints(int x, int y);
RcppExport SEXP _SparsePortfolioSelection_add_ints(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(add_ints(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SparsePortfolioSelection_add_ints", (DL_FUNC) &_SparsePortfolioSelection_add_ints, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SparsePortfolioSelection(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
