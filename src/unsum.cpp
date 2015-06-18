#include <Rcpp.h>

#include "dense_hybrid.h"

using namespace Rcpp;

void report(int result) {
  if(result) {
    switch(result) {
      case DH_FAIL_TIME_OUT: {
        Rcout << "  Quench timed out: Restarting..." << endl;
      } break;
      case DH_FAIL_PLATEAU: {
        Rcout << "  Energy plateaued: Restarting..." << endl;
      } break;
      case DH_FAIL_ROWCOL: {
        Rcout << "  Invalid minimum: Restarting..." << endl;
      } break;
    } // End switch
  }
}

// This function faces R. It creates a dense_hybrid object and launches
// the simulations, extracting the results and sending back to R.

// [[Rcpp::export]]
SEXP unsumcpp(Rcpp::NumericMatrix constraints,
             int target_ne,
             bool maxEdges,
             bool noReturn,
             long  mct_schedule,
             long  hot_time,
             double beta0,
             double betamax,
             double mu0,
             double cooling_rate,
             long max_time,
             double cgmax)
{
    int nn = constraints.nrow();
    dense_hybrid *dh; // everything happens in this object
    
    dh = new dense_hybrid(nn, target_ne, maxEdges, noReturn);
    
    // Convert the constraints into something dense_hybrid can read
    // Out is first column, in is second column
    // rowSum==Out
    // colSum==In
    NumericVector tv = constraints(_,0);
    std::vector<double> targetsOut = as<std::vector<double> >(tv);
    tv = constraints(_,1);
    std::vector<double> targetsIn = as<std::vector<double> >(tv);
    
    // Copy in the constraints to the targets
    dh->init_targets(targetsOut, targetsIn);
    
    // You'll need to do this a lot more than once
    int result=1;
    
    while(result>0) {
      result = dh->runjob(mct_schedule, hot_time, beta0, betamax, mu0,
      cooling_rate, max_time, cgmax);
      
      report(result); // If result>0 then report back.
      
      checkUserInterrupt(); // This should happen in dense_hybrid but to be sure
      
      Rcout << "Result=" << result << endl;
    }
        
    IntegerMatrix A(nn, nn);
    NumericMatrix W(nn, nn);
    // Copy data from the matrices in dh
    // Both are column major (I think)
    int ne=0; // count the edges
    for (int i=0; i<nn*nn; i++) {
      A[i] = dh->A[i];
      ne += A[i];
      W[i] = dh->W[i] * dh->scalew;
    }
    
    // The refined solution is in the activeW and active_edges arrays
    NumericMatrix AW(nn, nn);
    for (int i=0; i<ne; i++) {
      int k = dh->active_edges[i];
      AW[k] = dh->activeW[i] * dh->scalew;
    }
    
    delete dh; // cleanup

    return List::create(Named( "A" ) = A,
                        Named( "W" ) = W,
                        Named( "AW") = AW);
}

