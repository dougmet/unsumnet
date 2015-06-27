#include "rng_wrap.h"

#ifdef BUILD_FOR_R

rng_wrap::rng_wrap(int oneSeed) {
    seed(oneSeed);
}

rng_wrap::rng_wrap() {
    // Nothing to do
}

double rng_wrap::rand53() {
    return(R::runif(0,1));
}

int rng_wrap::randInt(int N) {
    return(int((N+1)*R::runif(0,1))); // I don't think this can return N+1
}

void rng_wrap::seed(int oneSeed) {
    Rcpp::stop("Have not implemented Cpp-side one-seed setting yet.");
}
//void rng_wrap::seed(const uint32 * seedArray, int N) {
//    Rcpp::stop("Have not implemented Cpp-side big-seed setting yet.");
//}

#else

rng_wrap::rng_wrap(int oneSeed) {
    mt.seed(oneSeed);
}

rng_wrap::rng_wrap() {
    // Nothing to do
}

double rng_wrap::rand53() {
    return(mt.rand53());
}

int rng_wrap::randInt(int N) {
    return(mt.randInt(N));
}

void rng_wrap::seed(int oneSeed) {
    mt.seed(oneSeed);
}
//void rng_wrap::seed(const unsigned int * seedArray, int N) {
//    mt.seed(seedArray, N);
//}

#endif

