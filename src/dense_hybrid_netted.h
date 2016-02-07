/*
 
 This file is part of dense_hybrid_netted.
 
 dense_hybrid_netted is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 2 of the License, or
 (at your option) any later version.
 
 dense_hybrid_netted is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with dense_hybrid.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#ifndef DENSE_HYBRID_NETTED_H // header protection
#define DENSE_HYBRID_NETTED_H

#define DH_SUCCESS 0
#define DH_FAIL_TIME_OUT 1
#define DH_FAIL_PLATEAU 2
#define DH_FAIL_ROWCOL 3

// This code is for dense graphs so we have the entire adjancy matrix

#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>
#include <cstring>
using namespace std;

#ifdef BUILD_FOR_R
#include <Rcpp.h>
#else
// Nothing for now
#endif

#define FTOL 1e-16
//#define ITMAX 100000
#define ITMAX 10000

#include "MersenneTwister.h"  // Random numbers
#include "move_class.h"       // MC move selection and control



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
    int *bankid;        // A unique identifier
    double *target_out; // Target for sum of out weights
    double *target_in;  // Target for sum of in weights
    double *sum_out;    // Current sum of out weights
    double *sum_in;     // Current sum of in weights
    
    // Netted constraints
    double *target_net_out; // Target for netted sum out
    double *target_net_in;  // Target for netted sum in
    double *target_ex_out; // Target for netted sum to the external bank
    double *sum_net_out;    // Current netted sum out
    double *sum_net_in;     // Current netted sum in
    double *cg_sum_net_out;
    double *cg_sum_net_in;
    long *reverse_edges;     // Gives the index, if it exists, of the reverse edge in the active edges array
    
    bool *on_out;
    bool *on_in;
    bool *on_net_ex;
    double *top,*tip;   // Speed up the scaled energy calculation by storing some powers
    double *top_net,*tip_net;   // Speed up the scaled energy calculation by storing some powers
    double *top_ex;
    bool goforquench;
    
    // For conjugate gradient
    double *activeW;
    double *cg_sum_in;
    double *cg_sum_out;

    
    double energy, energynet, energyex;
    int nnet, nex;      // How many nodes are controlling for this
    
    
    MTRand mt;          // Random number generator

    double SCALE_FAC;

    double maxw,minw, scalew;    // restrict the range of Ws (careful!)
    long nn, ncycles, vol;
    long target_ne, ne, ne0, nepair;     // Target number of edges and actual number, number at end of hot time, number of edge pairs
    
    double beta, mu;        // fields for energy and ne (beta is inverse temperature, beta=1/T)
    double cooling_rate;
    
    
    ////////// METHODS //////////
    dense_hybrid(int nn_in, int target_ne_in,
                 bool verbose_in,
                 bool inMAXEDGES,
                 bool inNORETURN); // constructor
	~dense_hybrid();	// destructor
    
    int runjob(long  mct_schedule,
               long  hot_time,
               double beta0, double betamax,
               double betanet0, double betanetmax,
               double betaex0, double betaexmax,
               double mu0, double mumax,
               double cooling_rate,
               long max_time,
               double cgmax);
    
    int read_input(const char *infilename);
    void create_arrays();
    void initialise_arrays();
    void init_targets(const std::vector<double> & target_out_i,
                      const std::vector<double> & target_in_i,
                      const std::vector<double> & target_net_out_i,
                      const std::vector<double> & target_net_in_i,
                      const std::vector<double> & target_ex_out_i);
	void reset_arrays();
    double total_energy(double *inW);
    double total_energy();
    double total_energy_netted();
    double total_energy_ex();
    bool rowcol_iterate();
    void link_active_reverse();

    int mc_sweep(move_class ** move, int Nmoves);

    
};

#endif // end header protection
