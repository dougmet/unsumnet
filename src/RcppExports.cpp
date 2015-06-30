// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// maxEntropyCpp
NumericMatrix maxEntropyCpp(NumericMatrix aw, NumericVector rs, NumericVector cs, double minError);
RcppExport SEXP unsumnet_maxEntropyCpp(SEXP awSEXP, SEXP rsSEXP, SEXP csSEXP, SEXP minErrorSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type aw(awSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rs(rsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cs(csSEXP);
    Rcpp::traits::input_parameter< double >::type minError(minErrorSEXP);
    __result = Rcpp::wrap(maxEntropyCpp(aw, rs, cs, minError));
    return __result;
END_RCPP
}
// unsumcpp
SEXP unsumcpp(Rcpp::NumericMatrix constraints, int target_ne, bool maxEdges, bool noReturn, long mct_schedule, long hot_time, double beta0, double betamax, double mu0, double cooling_rate, long max_time, double cgmax);
RcppExport SEXP unsumnet_unsumcpp(SEXP constraintsSEXP, SEXP target_neSEXP, SEXP maxEdgesSEXP, SEXP noReturnSEXP, SEXP mct_scheduleSEXP, SEXP hot_timeSEXP, SEXP beta0SEXP, SEXP betamaxSEXP, SEXP mu0SEXP, SEXP cooling_rateSEXP, SEXP max_timeSEXP, SEXP cgmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< int >::type target_ne(target_neSEXP);
    Rcpp::traits::input_parameter< bool >::type maxEdges(maxEdgesSEXP);
    Rcpp::traits::input_parameter< bool >::type noReturn(noReturnSEXP);
    Rcpp::traits::input_parameter< long >::type mct_schedule(mct_scheduleSEXP);
    Rcpp::traits::input_parameter< long >::type hot_time(hot_timeSEXP);
    Rcpp::traits::input_parameter< double >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< double >::type betamax(betamaxSEXP);
    Rcpp::traits::input_parameter< double >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< double >::type cooling_rate(cooling_rateSEXP);
    Rcpp::traits::input_parameter< long >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< double >::type cgmax(cgmaxSEXP);
    __result = Rcpp::wrap(unsumcpp(constraints, target_ne, maxEdges, noReturn, mct_schedule, hot_time, beta0, betamax, mu0, cooling_rate, max_time, cgmax));
    return __result;
END_RCPP
}
