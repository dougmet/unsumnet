#include <Rcpp.h>
#include "dense_hybrid.h"

// This function faces R. It creates a dense_hybrid object and launches
// the simulations, extracting the results and sending back to R.

// [[Rcpp::export]]
int unsumcpp(Rcpp::NumericMatrix constraints,
             long  mct_schedule,
             long  hot_time,
             double beta0,
             double betamax,
             double mu0,
             double cooling_rate,
             long max_time,
             double cgmax)
{
    
    dense_hybrid *dh; // everything happens in this object
    
    dh = new dense_hybrid(nn, target_ne);
    
    // Convert the constraints into something dense_hybrid can read
    std::vector<double> targetsOut = Rcpp::as<std::vector<double> >(contstraints[SOMETHING]);
    std::vector<double> targetsIn = Rcpp::as<std::vector<double> >(contstraints[SOMETHING]);
    
    // Copy in the constraints to the targets
    dh->init_targets(targetsOut, targetsIn);
    
    // You'll need to do this a lot more than once
    dh->runjob(ncycles, mct_schedule, hot_time, beta0, betamax, mu0,
               cooling_rate, max_time, cgmax);
    
    delete dh;

    return 0;
}