#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
int timesTwo(int x) {
   return x * 2;
}

int main(int argc, char *argv[])
{
    
    dense_hybrid *dh; // everything happens in this object
    
    dh = new dense_hybrid;
    gdh = dh; // global pointer
    
    dh->runjob(argv[1]);

    return 0;
}