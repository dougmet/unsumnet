// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// unsumcpp
int unsumcpp(int awhole, int ton, double of, double variables);
RcppExport SEXP banknet_unsumcpp(SEXP awholeSEXP, SEXP tonSEXP, SEXP ofSEXP, SEXP variablesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type awhole(awholeSEXP);
    Rcpp::traits::input_parameter< int >::type ton(tonSEXP);
    Rcpp::traits::input_parameter< double >::type of(ofSEXP);
    Rcpp::traits::input_parameter< double >::type variables(variablesSEXP);
    __result = Rcpp::wrap(unsumcpp(awhole, ton, of, variables));
    return __result;
END_RCPP
}