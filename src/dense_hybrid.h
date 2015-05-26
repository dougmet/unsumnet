

#ifndef DENSE_HYBRID_H // header protection
#define DENSE_HYBRID_H

// This code is for dense graphs so we have the entire adjancy matrix

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>
#include <cstring>
using namespace std; 
#include "MersenneTwister.h"

#define FTOL 1e-16
//#define ITMAX 100000
#define ITMAX 10000


#ifdef ENERGYBIAS
  #include "tmclass.h"
#endif 

//#include "minimise.cpp" // No conjugate gradient for now

#include "move_class.h"



void write_monitor(const char *filename, long t, double data)
{
    ofstream monfile;
    monfile.open(filename, ios::app);
    monfile << t << " " << data << endl;
    monfile.close();
}

class dense_hybrid
{
    
public:
	//////// VARIABLES //////////
    
    // Order N*N arrays
    double *W;          // Edge weights
    int *A;             // Adjancancy Matrix
    long *active_edges; // Which edges are in use? Useful for sparser matrices
    long *active_edges0;// For correlation function
    
    // Order N arrays
    double *target_out; // Target for sum of out weights
    double *target_in;  // Target for sum of in weights
    double *sum_out;    // Current sum of out weights
    double *sum_in;     // Current sum of in weights
    bool *on_out;
    bool *on_in;
    double *top,*tip;   // Speed up the scaled energy calculation by storing some powers
    bool goforquench;
    
    // For conjugate gradient
    double *gradW;
    double *activeW;
    double *cg_sum_in;
    double *cg_sum_out;
    double *cg_g, *cg_h;

    double energy;
    
    MTRand mt;          // Random number generator

    double SCALE_FAC;

    double maxw,minw, scalew;    // restrict the range of Ws (careful!)
    long nn, ncycles, vol;
    long target_ne, ne, ne0;     // Target number of edges and actual number, number at end of hot time
    
    double beta, mu;        // fields for energy and ne (beta is inverse temperature, beta=1/T)
    double cooling_rate;
    
    
    ////////// METHODS (AKA FUNCTIONS) //////////
    dense_hybrid(int nn_in, int target_ne_in); // constructor
    int runjob(int ncycles,
                long  mct_schedule,
                long  hot_time,
                double beta0,
                double betamax,
                double mu,
                double cooling_rate,
                long max_time,
                double cgmax, // when to attempt conjugate gradient
                double CG_TARGET);

    int read_input(const char *infilename);
    void create_arrays();
    void initialise_arrays();
	void reset_arrays();
    double total_energy(double *inW);
    double total_energy();
    bool rowcol_iterate();
    void gradient_energy(double *inW, double *xi);
    void conjugate_gradient();
    
};

#endif // end header protection