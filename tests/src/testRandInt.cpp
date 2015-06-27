#include <Rcpp.h>
using namespace Rcpp;

int randIntCpp(int N) {
  return(int((N+1)*R::runif(0,1))); // I don't think this can return N+1
}

// [[Rcpp::export]]
bool testRandInt(int N) {

  for(int i=0; i< N; i++) {
    if (randIntCpp(10)==11)
      stop("randInt has returned out of range value");
  }
  
  return(true);
}

/*** R
for (i in 1:20) testRandInt(2147483647L)
*/
