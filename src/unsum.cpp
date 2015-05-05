#include <Rcpp.h>
#include "dense_hybrid.h"
using namespace Rcpp;


// [[Rcpp::export]]
int unsumcpp(int awhole,
          int ton,
          double of,
          double variables)
{
    
    dense_hybrid *dh; // everything happens in this object
    
    dh = new dense_hybrid(nn, target_ne);
    gdh = dh; // global pointer
    
    dh->runjob(awhole, ton, of, variables);
    
    delete dh;

    return 0;
}