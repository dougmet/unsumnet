
#ifndef RNG_WRAP_H // header protection
#define RNG_WRAP_H

#ifdef BUILD_FOR_R
 #include <Rcpp.h>
#else
 #include "MersenneTwister.h"
#endif


class rng_wrap
{
    
public:
	//////// VARIABLES //////////
#ifdef BUILD_FOR_R
    // Nothing to do
#else
    MTRand mt;
#endif

    ////////// METHODS (AKA FUNCTIONS) //////////
    rng_wrap(int oneSeed);  // constructors
    rng_wrap();
    
    double rand53();
    int randInt(int N);
    
    void seed(int oneSeed);
   // void seed(const unsigned int * seedArray, int N); // Removed for now
    
};

#endif // end header protection