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
 along with dense_hybrid_netted.  If not, see <http://www.gnu.org/licenses/>.
 
 */

// This code is for dense graphs so we have the entire adjancy matrix

// In this netted version the external bank MUST be node 0

#include "dense_hybrid_netted.h"


inline double pow2(double x)
{
    return(x*x);
}


dense_hybrid_netted::dense_hybrid_netted(int nn_in, int target_ne_in,
                           bool verbose_in,
                           bool inMAXEDGES,
                           bool inNORETURN)
{
    ///////////////////////////
    //// CREATE THE ARRAYS ////
    ///////////////////////////
    
    verbose = verbose_in;
    MAXEDGES = inMAXEDGES;
    NORETURN = inNORETURN;
    
    nn = nn_in;
    target_ne = target_ne_in;
    
    // Set the size of the edge arrays
    if (nn<=32)
        max_ne = nn*nn;
    else if (2*target_ne < 1024)
        max_ne = 1024;
    else
        max_ne = 2*target_ne;
    
    ///////////////////////////
    //// CREATE THE ARRAYS ////
    ///////////////////////////
    
    // They're 1D to keep pointers simple
    // I've put them in COLUMN-MAJOR format as this is compatible with blas and lapack
    // libraries incase we want eigenvectors etc.
    // To get element A_ij address A[i + j*nn]
    A = new int [nn*nn];
    W = new double [nn*nn];
    active_edges = new long [target_ne*2];            // first ne elements give active edges
    active_edges0 = new long [target_ne*2];    // For correlation functions
    vol = nn*(nn-1);                            // total number of possible edges
    
    // These arrays store node properties
    bankid = new int [nn];
    target_out = new double [nn];
    target_in = new double [nn];
    sum_out = new double [nn];
    sum_in = new double [nn];
    // Netted properties
    target_net_out = new double [nn];
    target_net_in = new double [nn];
    target_ex_out = new double [nn];
    sum_net_out = new double [nn];
    sum_net_in = new double [nn];
    top_net = new double [nn];
    tip_net = new double [nn];
    top_ex = new double [nn];
    reverse_edges = new long [target_ne*2];
    
    on_out = new bool [nn];
    on_in = new bool [nn];
    on_net_ex = new bool [nn];
    top = new double [nn];
    tip = new double [nn];
    
    // Arrays for total energy (formally conjugate gradient)
    activeW= new double [target_ne*2];
    cg_sum_out = new double [nn];
    cg_sum_in = new double [nn];
    cg_sum_net_out = new double [nn];
    cg_sum_net_in = new double [nn];
    
    
}


dense_hybrid_netted::~dense_hybrid_netted()
{
    // try to free everything
    delete[] A;
    delete[] W;
    delete[] active_edges;
    delete[] active_edges0;
    delete[] bankid;
    delete[] target_out;
    delete[] target_in;
    delete[] sum_out;
    delete[] sum_in;
    delete[] target_net_out;
    delete[] target_net_in;
    delete[] target_ex_out;
    delete[] sum_net_out;
    delete[] sum_net_in;
    delete[] top_net;
    delete[] tip_net;
    delete[] top_ex;
    delete[] reverse_edges;
    delete[] on_out;
    delete[] on_in;
    delete[] on_net_ex;
    delete[] top;
    delete[] tip;
    delete[] activeW;
    delete[] cg_sum_out;
    delete[] cg_sum_in;
    delete[] cg_sum_net_out;
    delete[] cg_sum_net_in;
    
}



void dense_hybrid_netted::init_targets(const std::vector<double> & target_out_i,
                                       const std::vector<double> & target_in_i,
                                       const std::vector<double> & target_net_out_i,
                                       const std::vector<double> & target_net_in_i,
                                       const std::vector<double> & target_ex_out_i)
{
    
    
    maxw=0;
    int ncount=0, checkid;
    
    
    for (int i=0;i<nn;i++)
    {
        
        on_out[i] = on_in[i] = on_net_ex[i] = true;
        if (fabs(target_out_i[i])<=1e-12)
            on_out[i]=false;
        if (fabs(target_in_i[i])<=1e-12)
            on_in[i]=false;
        
        
        // Netted to external is positive or negative
        if (target_ex_out_i[i]>-1e-12 && target_ex_out_i[i]<1e-12)
        {
            on_net_ex[i]=false;
        }
        
        if (target_in_i[i]>maxw)
            maxw = target_in_i[i];
        
        if (target_out_i[i]>maxw)
            maxw = target_out_i[i];
        
    }
    
    
    scalew = maxw;
    maxw=2.0;
    minw=0.01/scalew;
    
    // Now scale the whole problem
    for (int i=0;i<nn;i++)
    {
        target_out[i] = target_out_i[i] / scalew;
        target_in[i] = target_in_i[i] / scalew;
        target_net_out[i] = target_net_out_i[i] / scalew;
        target_net_in[i] = target_net_in_i[i] / scalew;
        target_ex_out[i] = target_ex_out_i[i] / scalew;
    }

    
}


int dense_hybrid_netted::runjob(long  mct_schedule,
                         long  hot_time,
                         double beta0, double betamax,
                         double betanet0, double betanetmax,
                         double betaex0, double betaexmax,
                         double mu0, double mumax,
                         double cooling_rate,
                         long max_time,
                         double cgmax)
{
    
    HERE HERE HERE

    // This function is the main job controller
    
    int newi,newj,nrn,newk;                          // Edge moves
    int ir, mr;                                 // for choosing moves
    long mct, mct_schedule, hot_time, max_time;       // Monte Carlo time
	
	double oldenergy;
	int nsmallderiv=0;
	
    bool targetsread=false;
    bool netout;        // useful for remembering which way the netted flow is
	

    
//    infile >> SCALE_FAC;
    SCALE_FAC=1.0;
    
    beta=beta0;
    betaex=betaex0;
	betanet=betanet0;
    mu=mu0;
    
    
    ///////////////////////////
    ///// INITIALISATION //////
    ///////////////////////////
    

    
	reset_arrays(); // this function zeroes A, resets W and calculates energy

    
    ////////// MOVES ////////////
	// Setup the moves we'll use
    const int Nmoves = 7;
    int moves_perMC, sum_moves;
    move_class *move[Nmoves];
    
    // Initialise the moves
    // Arguments are: Name, Moves per sweep, initial step, minimum step, maximum step, target acceptance rate
    move[0] = new move_class("Edge insertion", target_ne, 0, 0, 0, 1);       // no step size needed here
    move[1] = new move_class("Edge deletion", target_ne, 0, 0, 0, 1);        // no step size needed here
	move[2] = new move_class("Edge tweak", target_ne*4, 1.0, 0.0000001, 10.0, 0.40);   // Aiming for 40%
    move[3] = new move_class("Edge out rewire", 0*target_ne, 0, 0, 0, 1);        // no step size needed here
    move[4] = new move_class("Edge in rewire", 0*target_ne, 0, 0, 0, 1);        // no step size needed here
    move[5] = new move_class("Edge pair insertion", target_ne, 0, 0, 0, 1);       // no step size needed here
    move[6] = new move_class("Edge pair deletion", target_ne, 0, 0, 0, 1);        // no step size needed here
    if (MAXEDGES)     // In the max edges run we always switch off insertions/deletions and edge moves
    {
        move[0]->NperMC=move[1]->NperMC=move[3]->NperMC=move[4]->NperMC=0;
        move[5]->NperMC=move[6]->NperMC=0;
    }
    
    
    // Clear some files
    // Make a function for this
//    monfile.open("energy.dat"); monfile.close();
//    monfile.open("energynet.dat"); monfile.close();
//    monfile.open("energyex.dat"); monfile.close();
//    monfile.open("beta.dat"); monfile.close();
//    monfile.open("step.dat"); monfile.close();
//    monfile.open("success.dat"); monfile.close();
//    monfile.open("ne.dat"); monfile.close();
//    monfile.open("edge_corr.dat"); monfile.close();
//    monfile.open("maxlocal.dat"); monfile.close();
    
    
    ////////////////////////////////////
    //////////// MAIN LOOP /////////////
    ////////////////////////////////////
    
    for (mct=0;mct<max_time;mct++)
    {
        if(!MAXEDGES)    // In the max edges run we always switch off insertions/deletions and edge moves
        {
            if (!goforquench) // restore insertions/deletions
                move[0]->NperMC=move[1]->NperMC=move[5]->NperMC=move[6]->NperMC=target_ne;
        }
        else
        {
            goforquench=true;
        }
        
        ////////////////////////////
        ////    MAIN MC LOOP    ////
        ////////////////////////////
        

        
         ///////////////////////////////////////////////////////////
                        //// MEASUREMENTS ////
         ///////////////////////////////////////////////////////////
        
        if ((mct%mct_schedule)==0)
        {
            // Change the temperature in the annealing scheme
            if (mct>hot_time)
            {
                // Cooling strategy
                beta *= cooling_rate;
                if (mu<mumax)
                    mu *= cooling_rate;
                if (betanet < betanetmax)
                    betanet *= cooling_rate;
                if (betaex < betaexmax)
                    betaex *= cooling_rate;

            }
            
            // Update the step sizes based on the acceptance rate
            for (mr=0;mr<Nmoves;mr++)
				move[mr]->update_step();
            
			if ((mct%(100*mct_schedule))==0)
            {
                
#ifdef BUILD_FOR_R
                if(verbose) {

                    Rcout <<  "  t=" << mct << ", beta=" << beta << ", step=" <<
                      move[2]->step << ":" << move[2]->success_rate << ", E=" <<
                      energy/nn << ", Enet=" << energynet/nnet << ", Eex="<<
                      energyex/nex <<", ne=" << ne << ":" << move[0]->success_rate <<
                      ":" << move[1]->success_rate << endl;
                      //" rw" << move[3]->success_rate << ":" << move[4]->success_rate << endl;
                }
#endif
            

                energy = total_energy();
                energynet = total_energy_netted();
                energyex = total_energy_ex();

                if (mct>1000)
                {
                    write_monitor("energy.dat", mct, energy/nn);
                    write_monitor("energynet.dat", mct, energynet/nnet);
                    write_monitor("energyex.dat", mct, energyex/nex);
                    write_monitor("beta.dat", mct, beta);
                    //write_monitor("step.dat", mct, move[2]->step);
                    //write_monitor("success.dat", mct, move[2]->success_rate);
                    write_monitor("ne.dat", mct, (double) ne);
                }
				


				
            }

            
            ///////////////////////////////////////////////////////////
                /////////////////  QUENCH CONTROL ///////////////
            ///////////////////////////////////////////////////////////
            
            if (mct<=hot_time)
            {
                // Save a list of who is switched on before we start cooling
                ne0=ne;
                for (i=0;i<ne;i++)
                    active_edges0[i]=active_edges[i];
				
				if (mct==hot_time)
				{
#ifdef BUILD_FOR_R
                    if(verbose)
                        Rcout << "End of hot_time" << endl;
#endif
					energy = total_energy();
				}
                
            }
            else
            {
              
   
                for (k=0, i=0;i<ne0;i++)
                    k += A[active_edges0[i]]; // sum all edges who are still on
                write_monitor("edge_corr.dat", mct, ((double) k)/ne0);

                // Do you want to abort this quench?
				if ((mct%(20*mct_schedule))==0)
				{
					if (fabs(1-oldenergy/energy)< 0.001)
						nsmallderiv++;
					else
						nsmallderiv=0;
					oldenergy=energy;
					
					if (nsmallderiv>10)
					{
                        // Plateau
                        return(DH_FAIL_PLATEAU);
					}
				}
				
                // Row-Column iteration solution finder
                if (mct%(50*mct_schedule)==0)
                {

					
                    if (energy/nn<5e-3)
                    if (move[0]->success_rate < 1e-3 && !goforquench)
                    {
                        if (!rowcol_iterate())
                        {
                            // Invalid structure
                            return(DH_FAIL_ROWCOL);
                        }
                        else
                        {
                            if (ne==target_ne)
                            {
#ifdef BUILD_FOR_R
                                if(verbose)
                                    Rcout << "Switching off insertions/deletions." << endl;
#endif
                                goforquench=true;
                                // Switch off insertions/deletions, we'll stick with this config
                                move[0]->NperMC=move[1]->NperMC=0;
                                move[5]->NperMC=move[6]->NperMC=0;
                            }
                        }
                    }

                    
                    // Kill edge rewiring when that stops
                    if (move[3]->success_rate < 1e-3 && move[3]->NperMC > 0.1)
                    {
                        if (rowcol_iterate())
                        {
#ifdef BUILD_FOR_R
                            if(verbose)
                                Rcout << "Freezing edge swaps." << endl;
#endif
                            move[3]->NperMC=move[4]->NperMC=0;
                        }
                    }
                    
                }
                
                
                ///////////////////////////////////////////////////////////
                    /////////////////  FINAL ANSWER ///////////////
                ///////////////////////////////////////////////////////////
                
                // Conjugate gradient
                if (mct%(50*mct_schedule)==0)
                if (energy<cgmax*nn && rowcol_iterate())
                {
                    // CLose enough
                    return(DH_SUCCESS);
                }

				
            }
            
#ifdef BUILD_FOR_R
            // Check in with R
            Rcpp::checkUserInterrupt();
#endif

        }
        
    } // end mct loop
	
    // Hit the buffers without success
    return (DH_FAIL_TIME_OUT);
}

//////////////////////////////////
//      Monte Carlo Sweep       //
//////////////////////////////////

int dense_hybrid_netted::mc_sweep(move_class ** move, int Nmoves)
{
    
    int moves_perMC, sum_moves;
    long i,j,k,kr,iact, irev;                           // general integers
    int newi,newj,nrn,newk;             // Edge moves
    int ir, mr;                         // for choosing moves
    double deltaE, deltaM, w, dw;
    double r, prob;
    double localmaxw;   // Working numbers
    double deltaE, deltaM, w, dw, r, prob;
    double tenergy, temax, localmaxw, wnet;    // Working numbers

    
    // This tells us how many moves per MC sweep there are in total
    for (moves_perMC=0, j=0;j<Nmoves;j++)
        moves_perMC += move[j]->NperMC;
    
    //// BEGIN MONTE CARLO SWEEP ////
    for (int im=0;im<moves_perMC;im++)
    {
        // Choose a move (with replacement)
        ir = mt.randInt(moves_perMC-1);     // pick a number from the total number
        for (sum_moves=0, mr=0; mr<Nmoves; mr++)
        {
            sum_moves+=move[mr]->NperMC;
            if (ir < sum_moves)
                break;
        }
        // when you get here mr corresponds to the move index
        
        
        switch (mr)
        {
                
            case 0:
            {
                ////// EDGE INSERTION //////
                
                // Choose random edge
                i = mt.randInt(nn-1); // from i
                do {
                    j = mt.randInt(nn-1); // to j
                } while (j==i); // no self loops
                
                k = i + j*nn;   // position in the adjacency and weights matrix
                kr = j + i*nn;  // return edge
                
                
                if (on_out[i] && on_in[j]) // This stops edges to nodes with a zero in/out sum
                {
                    if (A[k] == 0)
                    {
                        // Calculate the change at the out & in nodes
                        w = W[k];
                        
                        //////////////////////// GROSS SUMS //////////////////////////
                        deltaE = w*(w - 2*(target_out[i] - sum_out[i]))/pow2(top[i]);    // Out node
                        deltaE += w*(w - 2*(target_in[j] - sum_in[j]))/pow2(tip[j]);     // In node
                        
                        //////////////////////// NETTED SUMS //////////////////////////
                        if (A[kr])
                        {
                            // First remove the previous netted position
                            wnet = w - W[kr];   // What the netted position will become
                            if (wnet>0)
                            {
                                netout = true; // net direction is i to j
                                // flow changes direction, remove old flow
                                deltaEnet = W[kr]*(W[kr] + 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]);
                                deltaEnet += W[kr]*(W[kr] + 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]);
                                // now the new net is applied to the new direction
                                deltaEnet += wnet*(wnet - 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]); // outbound from i
                                deltaEnet += wnet*(wnet - 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]); // inbound to j
                            }
                            else
                            {
                                netout = false; // net direction is j to i
                                // flow retains its direction, subtract W[k] from its value
                                deltaEnet = W[k]*(W[k] + 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]); // inbound to i
                                deltaEnet += W[k]*(W[k] + 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]); // outbound from j
                            }
                            
                        }
                        else
                        {
                            // Flow is guaranteed from i to j as it's the only edge
                            wnet = w;
                            netout = true;
                            
                            deltaEnet = wnet*(wnet - 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]); // outbound from i
                            deltaEnet += wnet*(wnet - 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]); // inbound to j
                        }
                        
                        //////////////////////// EXTERNAL NET //////////////////////////
                        deltaEex=0;
                        if (i==0) // i is the external bank
                        {
                            if (on_net_ex[j])
                                deltaEex=w*(w + 2*(target_ex_out[j] - A[kr]*W[kr]))/pow2(top_ex[j]);
                        }
                        if (j==0) // j is the external bank
                        {
                            if (on_net_ex[i])
                                deltaEex=w*(w - 2*(target_ex_out[i] + A[kr]*W[kr]))/pow2(top_ex[i]);
                        }
                        
                        
                        //////////////////////// SPARSITY //////////////////////////
                        deltaM = 1 - 2*(target_ne - ne);
                        
                        // Accept or reject
                        prob = exp(-beta*deltaE - mu*deltaM - deltaEnet*betanet - deltaEex*betaex) * vol / ((double) ne + 1.0);
                        
                        
                        if (mt() < prob)
                        {
                            // Success!
                            A[k] = 1;                // switch on the edge
                            
                            active_edges[ne] = k;    // Add it to the end of the active list
                            ne ++;                          // Increment number of edges
                            
                            sum_out[i] += w; // update the totals
                            sum_in[j] += w;
                            energy += deltaE;
                            
                            // netted positions
                            energynet += deltaEnet;
                            
                            if (A[kr])
                            {
                                nepair ++; // count the return pairs
                                
                                if (netout)
                                {
                                    sum_net_out[j] -= W[kr];
                                    sum_net_in[i] -= W[kr];
                                    sum_net_out[i] += wnet;
                                    sum_net_in[j] += wnet;
                                }
                                else
                                {
                                    sum_net_in[i] -= w; // w=W[k]
                                    sum_net_out[j] -= w;
                                }
                                
                            }
                            else
                            {
                                sum_net_out[i] += wnet;
                                sum_net_in[j] += wnet;
                            }
                            
                            // External bank
                            energyex += deltaEex;
                            
                            move[mr]->update_attempts(1); // move sucessful
                        }
                        else
                        {
                            move[mr]->update_attempts(0); // move unsucessful
                        }
                        
                        
                    }
                    else
                    {
                        move[mr]->update_attempts(0); // move unsucessful
                    }
                }
                
            }
                break;
                
            case 1:
            {
                ///// EDGE DELETION //////
                
                if (ne>0)
                {
                    // Choose random edge from the active list
                    iact = mt.randInt(ne-1);    // iact is the position in the active list
                    k = active_edges[iact];     // k is the position in the adjacency and weights matrices
                    
                    // Resolve this into the from (i) and to edge (j)
                    j = k/nn;
                    i = k - j*nn;
                    
                    kr = j + i*nn;  // return edge
                    
                    // What's the edge weight for the chosen edge
                    w = W[k];
                    
                    ///////////////////////// GROSS SUMS //////////////////////////
                    deltaE = w*(w + 2*(target_out[i] - sum_out[i]))/pow2(top[i]);    // Out node
                    deltaE += w*(w + 2*(target_in[j] - sum_in[j]))/pow2(tip[j]);     // In node
                    
                    //////////////////////// NETTED SUMS //////////////////////////
                    if (A[kr])
                    {
                        wnet = W[k] - W[kr];
                        if (wnet>0)
                        {
                            // It flowed from i to j so we remove |wnet| from out[i] and in[j] so the direction is changing
                            netout = true;
                            deltaEnet = wnet*(wnet + 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]);
                            deltaEnet += wnet*(wnet + 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]);
                            // and then we add W[kr] to out[j] and in[i]
                            deltaEnet += W[kr]*(W[kr] - 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]);
                            deltaEnet += W[kr]*(W[kr] - 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]);
                        }
                        else
                        {
                            // It flowed from j to i so direction stays the same
                            // simply add W[k] to the totals
                            netout = false;
                            deltaEnet = W[k]*(W[k] - 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]);
                            deltaEnet += W[k]*(W[k] - 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]);
                        }
                    }
                    else
                    {
                        // Only need to worry about one edge
                        wnet = W[k];
                        netout=true;
                        deltaEnet = wnet*(wnet + 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]);
                        deltaEnet += wnet*(wnet + 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]);
                    }
                    // End netted
                    
                    //////////////////////// EXTERNAL NET //////////////////////////
                    deltaEex=0;
                    if (i==0) // i is the external bank
                    {
                        if (on_net_ex[j])
                            deltaEex= -w*(w + 2*(target_ex_out[j] - A[kr]*W[kr]))/pow2(top_ex[j]);
                    }
                    if (j==0) // j is the external bank
                    {
                        if (on_net_ex[i])
                            deltaEex= -w*(w - 2*(target_ex_out[i] + A[kr]*W[kr]))/pow2(top_ex[i]);
                    }
                    
                    
                    
                    //////////////////////// SPARSITY //////////////////////////
                    deltaM = 1 + 2*(target_ne - ne);
                    
                    
                    // Accept or reject
                    prob = exp(-beta*deltaE - mu*deltaM - deltaEnet*betanet - deltaEex*betaex) * ((double) ne) / vol;
                    
                    
                    if (mt() < prob)
                    {
                        // Success!
                        A[k] = 0;                               // Switch off the edge
                        active_edges[iact]=active_edges[ne-1];  // Replace the entry in the active list with the one at the end
                        ne --;                                  // Decrement the number of edges
                        
                        sum_out[i] -= w; // update the totals
                        sum_in[j] -= w;
                        
                        energy += deltaE;
                        
                        // netted positions
                        if (A[kr])
                        {
                            nepair --; // count the return pairs
                            
                            if (netout)
                            {
                                // It flowed from i to j so we remove |wnet| from out[i] and in[j] so the direction is changing
                                sum_net_out[i] -= wnet;
                                sum_net_in[j] -= wnet;
                                sum_net_out[j] += W[kr];
                                sum_net_in[i] += W[kr];
                            }
                            else
                            {
                                // It flowed from j to i so direction stays the same
                                // simply add W[k] to the totals
                                sum_net_in[i] += W[k];
                                sum_net_out[j] += W[k];
                            }
                        }
                        else
                        {
                            // Only need to worry about one edge
                            sum_net_out[i] -= wnet;
                            sum_net_in[j] -= wnet;
                        }
                        energynet += deltaEnet;
                        // end netter positions
                        
                        // External bank netting
                        energyex += deltaEex;
                        
                        move[mr]->update_attempts(1); // Move sucessful
                    }
                    else
                    {
                        move[mr]->update_attempts(0); // Move unsuccessful
                    }
                    
                }
                
                
            }
                break;
                
            case 2:
            {
                ///// EDGE TWEAK //////
                
                if (ne>0)
                {
                    // Choose random edge from the active list
                    iact = mt.randInt(ne-1);    // iact is the position in the active list
                    k = active_edges[iact];     // k is the position in the adjacency and weights matrices
                    
                    // Resolve this into the from (i) and to edge (j)
                    j = k/nn;
                    i = k - j*nn;
                    kr = j + i*nn;  // return edge
                    
                    r = move[mr]->step*(1.0 - 2*mt());    // evenly distributed on +/- stepsize
                    
                    if (target_out[i]>target_in[j])
                        localmaxw=1.2*target_in[j];
                    else
                        localmaxw=1.2*target_out[i];
                    
                    if (localmaxw<-1e-15)
                        localmaxw = maxw;
                    
                    dw = r*localmaxw;
                    
                    
                    if ((W[k] + dw > minw) && (W[k] + dw < localmaxw))
                    {
                        deltaE = dw*(dw - 2*(target_out[i] - sum_out[i]))/pow2(top[i]);    // Out node
                        deltaE += dw*(dw - 2*(target_in[j] - sum_in[j]))/pow2(tip[j]);     // In node
                        
                        // Calculate change to the netted positions
                        if (A[kr])
                        {
                            wnet = W[k] - W[kr];
                            if (wnet>0)
                            {
                                if (wnet+dw>0)
                                {
                                    // same direction
                                    deltaEnet = dw*(dw - 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]);     // Out node
                                    deltaEnet += dw*(dw - 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]);      // In node
                                }
                                else
                                {
                                    // flow changes direction, remove wnet and add new flows [ (wnet+dw) <0 ]
                                    deltaEnet = wnet*(wnet + 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]);     // Out node
                                    deltaEnet += wnet*(wnet + 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]);      // In node
                                    deltaEnet += (wnet+dw)*((wnet+dw) + 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]);     // Out node
                                    deltaEnet += (wnet+dw)*((wnet+dw) + 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]);      // In node
                                }
                            }
                            else
                            {
                                if (wnet + dw<0)
                                {
                                    // subtract dw from opposite directions
                                    deltaEnet = dw*(dw + 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]);     // Out node
                                    deltaEnet += dw*(dw + 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]);      // In node
                                }
                                else
                                {
                                    // flow changes direction, remove wnet(-ve) and add new flows
                                    deltaEnet = wnet*(wnet - 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]);     // Out node
                                    deltaEnet += wnet*(wnet - 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]);      // In node
                                    deltaEnet += (wnet+dw)*((wnet+dw) - 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]);     // Out node
                                    deltaEnet += (wnet+dw)*((wnet+dw) - 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]);      // In node
                                }
                                
                            }
                        }
                        else
                        {
                            // Flow stays same direction by default
                            deltaEnet = dw*(dw - 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]);     // Out node
                            deltaEnet += dw*(dw - 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]);      // In node
                        }
                        
                        
                        //////////////////////// EXTERNAL NET //////////////////////////
                        deltaEex=0;
                        if (i==0) // i is the external bank
                        {
                            if (on_net_ex[j])
                                deltaEex= dw*(dw + 2*(target_ex_out[j] + (W[k] - A[kr]*W[kr])))/pow2(top_ex[j]);
                            
                        }
                        if (j==0) // j is the external bank
                        {
                            if (on_net_ex[i])
                                deltaEex= dw*(dw - 2*(target_ex_out[i] - (W[k] - A[kr]*W[kr])))/pow2(top_ex[i]);
                        }
                        
                        
                        
                        // Accept or reject
                        prob = exp(-(deltaE*beta + deltaEnet*betanet + deltaEex*betaex));
                        
                        if (mt() < prob)
                        {
                            W[k] += dw;
                            sum_out[i] += dw;
                            sum_in[j] += dw;
                            energy += deltaE;
                            // Now the netted positions
                            if (A[kr])
                            {
                                if (wnet>0)
                                {
                                    if (wnet+dw>0)
                                    {   // stay same i to j
                                        sum_net_out[i] += dw;
                                        sum_net_in[j] += dw;
                                    }
                                    else
                                    {
                                        // flow changes direction
                                        sum_net_out[i] -= wnet;
                                        sum_net_in[j] -= wnet;
                                        sum_net_out[j] -= wnet+dw;
                                        sum_net_in[i] -= wnet + dw;
                                    }
                                }
                                else
                                {
                                    if (wnet + dw<0)
                                    {   // stay same j to i
                                        sum_net_out[j] -= dw;
                                        sum_net_in[i] -= dw;
                                    }
                                    else
                                    {
                                        // flow changes direction (wnet is -ve)
                                        sum_net_out[j] += wnet;
                                        sum_net_in[i] += wnet;
                                        sum_net_out[i] += wnet+dw;
                                        sum_net_in[j] += wnet+dw;
                                    }
                                    
                                }
                            }
                            else
                            {
                                sum_net_out[i] += dw;
                                sum_net_in[j] += dw;
                            }
                            energynet += deltaEnet;
                            // end netted
                            
                            // External bank netting
                            energyex += deltaEex;
                            
                            move[mr]->update_attempts(1); // move sucessful
                            
                        }
                        else
                        {
                            move[mr]->update_attempts(0); // move unsucessful
                        }
                    }
                    else
                    {
                        move[mr]->update_attempts(0);
                    }
                    
                }
                
                
            }
                break;
                
            case 3:
            {
                ///////// REWIRE OUT EDGE /////////
                cout << "NOT CODED YET" << endl; exit(1);
                if (ne>0)
                {
                    // Choose random edge from the active list
                    iact = mt.randInt(ne-1);    // iact is the position in the active list
                    k = active_edges[iact];     // k is the position in the adjacency and weights matrices
                    
                    // Resolve this into the from (i) and to edge (j)
                    j = k/nn;
                    i = k - j*nn;
                    
                    // What's the edge weight for the chosen edge
                    w = W[k];
                    
                    // Now choose another node to wire to
                    nrn=0;
                    do {
                        newj = mt.randInt(nn-1);
                        nrn++;
                        // Must not be on, self looping or deactivated
                    } while ((newj==i || newj==j || (!on_in[newj]) || A[i + newj*nn]) && (nrn<nn));
                    
                    if (target_out[i]>target_in[newj])
                        localmaxw=1.2*target_in[newj];
                    else
                        localmaxw=1.2*target_out[i];
                    
                    if (localmaxw<-1e-15)
                        localmaxw = maxw;
                    
                    if (nrn<nn && w<localmaxw)
                    {
                        
                        deltaE = w*(w - 2*(target_in[newj] - sum_in[newj]))/pow2(tip[newj]);   // New in node
                        deltaE += w*(w + 2*(target_in[j] - sum_in[j]))/pow2(tip[j]);           // Old in node
                        
                        //////// netted positions
                        
                        //////// end netted
                        
                        // Accept or reject
                        prob = exp(-deltaE*beta);
                        
                        
                        if (mt() < prob)
                        {
                            sum_in[j] -= w;
                            sum_in[newj] += w;
                            energy += deltaE;
                            move[mr]->update_attempts(1); // move sucessful
                            newk = i + newj*nn;
                            // Swap the edge weights around
                            W[k] = W[newk];
                            W[newk] = w;
                            A[k] = 0;
                            A[newk]=1;
                            // Update active list
                            active_edges[iact] = newk;
                            
                        }
                        else
                        {
                            move[mr]->update_attempts(0); // move unsucessful
                        }
                    }
                    else
                    {
                        move[mr]->update_attempts(0); // move unsucessful
                    }
                    
                }
                
            }
                break;
                
            case 4:
            {
                ///////// REWIRE IN EDGE /////////
                cout << "NOT CODED YET" << endl; exit(1);
                if (ne>0)
                {
                    // Choose random edge from the active list
                    iact = mt.randInt(ne-1);    // iact is the position in the active list
                    k = active_edges[iact];     // k is the position in the adjacency and weights matrices
                    
                    // Resolve this into the from (i) and to edge (j)
                    j = k/nn;
                    i = k - j*nn;
                    
                    // What's the edge weight for the chosen edge
                    w = W[k];
                    
                    // Now choose another node to wire to
                    nrn=0;
                    do {
                        newi = mt.randInt(nn-1);
                        nrn++;
                        // Must not be on, self looping or deactivated
                    } while ((newi==j || newi==i || (!on_out[newi]) || A[newi + j*nn]) && (nrn<nn));
                    
                    if (target_out[newi]>target_in[j])
                        localmaxw=1.2*target_in[j];
                    else
                        localmaxw=1.2*target_out[newi];
                    
                    if (localmaxw<-1e-15)
                        localmaxw = maxw;
                    
                    if (nrn<nn && w<localmaxw)
                    {
                        deltaE = w*(w + 2*(target_out[i] - sum_out[i]))/pow2(top[i]);          // Old out node
                        deltaE += w*(w - 2*(target_out[newi] - sum_out[newi]))/pow2(top[newi]);// New out node
                        
                        // Accept or reject
                        prob = exp(-deltaE*beta);
                        
                        
                        if (mt() < prob)
                        {
                            sum_out[i] -= w;
                            sum_out[newi] += w;
                            energy += deltaE;
                            move[mr]->update_attempts(1); // move sucessful
                            newk = newi + j*nn;
                            // Swap the edge weights around
                            W[k] = W[newk];
                            W[newk] = w;
                            A[k] = 0;
                            A[newk]=1;
                            // Update active list
                            active_edges[iact] = newk;
                            
                        }
                        else
                        {
                            move[mr]->update_attempts(0); // move unsucessful
                        }
                    }
                    else
                    {
                        move[mr]->update_attempts(0); // move unsucessful
                    }
                    
                }
                
            }
                break;
                
                
            case 5:
            {
                ////// EDGE PAIR INSERTION //////
                
                // Choose random edge
                i = mt.randInt(nn-1); // from i
                do {
                    j = mt.randInt(nn-1); // to j
                } while (j==i); // no self loops
                
                k = i + j*nn;   // position in the adjacency and weights matrix
                kr = j + i*nn;  // return edge
                
                
                if (on_out[i] && on_in[j] && on_out[j] && on_in[i]) // This stops edges to nodes with a zero in/out sum
                {
                    if (A[k] == 0 && A[kr]==0)
                    {
                        //////////////////////// GROSS SUMS //////////////////////////
                        w = W[k];   // i -> j
                        deltaE = w*(w - 2*(target_out[i] - sum_out[i]))/pow2(top[i]);    // Out node
                        deltaE += w*(w - 2*(target_in[j] - sum_in[j]))/pow2(tip[j]);     // In node
                        w = W[kr];   // j -> i
                        deltaE += w*(w - 2*(target_out[j] - sum_out[j]))/pow2(top[j]);    // Out node
                        deltaE += w*(w - 2*(target_in[i] - sum_in[i]))/pow2(tip[i]);     // In node
                        
                        
                        //////////////////////// NETTED SUMS //////////////////////////
                        
                        wnet = W[k] - W[kr];   // What the netted position will become i->j
                        if (wnet>0)
                        {
                            netout = true; // net direction is i to j
                            deltaEnet = wnet*(wnet - 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]); // outbound from i
                            deltaEnet += wnet*(wnet - 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]); // inbound to j
                        }
                        else
                        {
                            netout = false; // net direction is j to i
                            deltaEnet = wnet*(wnet + 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]); // outbound from j
                            deltaEnet += wnet*(wnet + 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]); // inbound to i
                        }
                        
                        
                        //////////////////////// EXTERNAL NET //////////////////////////
                        deltaEex=0;
                        if (i==0) // i is the external bank
                        {
                            if (on_net_ex[j])
                                deltaEex=wnet*(wnet + 2*target_ex_out[j])/pow2(top_ex[j]);
                        }
                        if (j==0) // j is the external bank
                        {
                            if (on_net_ex[i])
                                deltaEex=wnet*(wnet - 2*target_ex_out[i])/pow2(top_ex[i]);
                        }
                        
                        //////////////////////// SPARSITY //////////////////////////
                        deltaM = 4*(1 + ne - target_ne);
                        
                        // Accept or reject
                        prob = exp(-beta*deltaE - mu*deltaM - deltaEnet*betanet - deltaEex*betaex) * vol / ((double) ne + 1.0);
                        
                        
                        if (mt() < prob)
                        {
                            // Success!
                            A[k] = A[kr] = 1;                // switch on the edge
                            
                            
                            active_edges[ne] = k;           // Add it to the end of the active list
                            ne ++;                          // Increment number of edges
                            active_edges[ne] = kr;          // Add it to the end of the active list
                            ne ++;                          // Increment number of edges
                            nepair ++;                      // count the return pairs
                            
                            
                            sum_out[i] += W[k]; // update the totals
                            sum_in[j] += W[k];
                            sum_out[j] += W[kr]; // update the totals
                            sum_in[i] += W[kr];
                            energy += deltaE;
                            
                            // netted positions
                            energynet += deltaEnet;
                            
                            if (netout)
                            {
                                sum_net_out[i] += wnet;
                                sum_net_in[j] += wnet;
                            }
                            else
                            {
                                sum_net_out[j] -= wnet;
                                sum_net_in[i] -= wnet;
                            }
                            
                            // External bank
                            energyex += deltaEex;
                            
                            move[mr]->update_attempts(1); // move sucessful
                        }
                        else
                        {
                            move[mr]->update_attempts(0); // move unsucessful
                        }
                        
                        
                    }
                    else
                    {
                        move[mr]->update_attempts(0); // move unsucessful
                    }
                }
                
            }
                break;
                
                
                
            case 6:
            {
                ///// EDGE PAIR DELETION //////
                
                if (ne>0)
                {
                    // Choose random edge from the active list
                    iact = mt.randInt(ne-1);    // iact is the position in the active list
                    k = active_edges[iact];     // k is the position in the adjacency and weights matrices
                    
                    // Resolve this into the from (i) and to edge (j)
                    j = k/nn;
                    i = k - j*nn;
                    
                    kr = j + i*nn;  // return edge
                    
                    if (A[kr])
                    {
                        
                        ///////////////////////// GROSS SUMS //////////////////////////
                        w=W[k];
                        deltaE = w*(w + 2*(target_out[i] - sum_out[i]))/pow2(top[i]);    // Out node
                        deltaE += w*(w + 2*(target_in[j] - sum_in[j]))/pow2(tip[j]);     // In node
                        w=W[kr];
                        deltaE += w*(w + 2*(target_out[j] - sum_out[j]))/pow2(top[j]);    // Out node
                        deltaE += w*(w + 2*(target_in[i] - sum_in[i]))/pow2(tip[i]);     // In node
                        //////////////////////// NETTED SUMS //////////////////////////
                        
                        wnet = W[k] - W[kr];
                        if (wnet>0)
                        {
                            // It flowed from i to j
                            netout = true;
                            deltaEnet = wnet*(wnet + 2*(target_net_out[i] - sum_net_out[i]))/pow2(top_net[i]);
                            deltaEnet += wnet*(wnet + 2*(target_net_in[j] - sum_net_in[j]))/pow2(tip_net[j]);
                        }
                        else
                        {
                            // It flowed from j to i
                            netout = false;
                            deltaEnet = wnet*(wnet - 2*(target_net_out[j] - sum_net_out[j]))/pow2(top_net[j]);
                            deltaEnet += wnet*(wnet - 2*(target_net_in[i] - sum_net_in[i]))/pow2(tip_net[i]);
                        }
                        // End netted
                        
                        //////////////////////// EXTERNAL NET //////////////////////////
                        deltaEex=0;
                        if (i==0) // i is the external bank
                        {
                            if (on_net_ex[j])
                                deltaEex= -wnet*(wnet + 2*target_ex_out[j])/pow2(top_ex[j]);
                        }
                        if (j==0) // j is the external bank
                        {
                            if (on_net_ex[i])
                                deltaEex= -wnet*(wnet - 2*target_ex_out[i])/pow2(top_ex[i]);
                        }
                        
                        
                        
                        //////////////////////// SPARSITY //////////////////////////
                        deltaM = 4*(1 - ne + target_ne);
                        
                        // Accept or reject
                        prob = exp(-beta*deltaE - mu*deltaM - deltaEnet*betanet - deltaEex*betaex) * ((double) ne) / vol;
                        
                        
                        if (mt() < prob)
                        {
                            // Success!
                            A[k] = 0;                       // Switch off the first edge
                            active_edges[iact]=active_edges[ne-1];  // Replace the entry in the active list with the one at the end
                            ne --;                                  // Decrement the number of edges
                            
                            if (iact<ne-1 && active_edges[iact+1]==kr)
                                irev = iact+1;
                            else if (iact>0 && active_edges[iact-1]==kr)
                                irev = iact-1;
                            else
                            {
                                irev=0; // A horrible search ensues
                                while (active_edges[irev] != kr)
                                {
                                    irev ++;
                                }
                            }
                            
                            A[kr] = 0;                      // Switch off the second edge
                            active_edges[irev]=active_edges[ne-1];  // Replace the entry in the active list with the one at the end
                            ne --;                                  // Decrement the number of edges
                            
                            sum_out[i] -= W[k]; // update the totals
                            sum_in[j] -= W[k];
                            sum_out[j] -= W[kr];
                            sum_in[i] -= W[kr];
                            
                            energy += deltaE;
                            
                            // netted positions
                            
                            nepair --; // count the return pairs
                            
                            if (netout)
                            {
                                // It flowed from i to j
                                sum_net_out[i] -= wnet;
                                sum_net_in[j] -= wnet;
                            }
                            else
                            {
                                // It flowed from j to i
                                sum_net_in[i] += wnet;
                                sum_net_out[j] += wnet;
                            }
                            
                            energynet += deltaEnet;
                            // end netter positions
                            
                            // External bank netting
                            energyex += deltaEex;
                            
                            move[mr]->update_attempts(1); // Move sucessful
                        }
                        else
                        {
                            move[mr]->update_attempts(0); // Move unsuccessful
                        }
                        
                    } // A[kr]
                } // ne>0
                
            }
                break;
                
                
        } // END SWITCH
        
        
    }
    //// END OF MONTE CARLO SWEEP ////
    
    
}

//////////////////////////////////
// Conjugate gradient functions
//////////////////////////////////

void dense_hybrid_netted::reset_arrays()
{
    goforquench=false;
    
	// Initialisation of edges
    for (int i=0;i<nn*nn;i++)
    {
        A[i]=0;         // All switched OFF
        W[i]=minw;   // All weights equal (could set to anything)
    }
    ne=0;               // No edges at the start
    nepair=0;           // No edge pairs at the start
    
    // Initialise the sums of edges
    for (int i=0;i<nn;i++)
    {
        sum_out[i] = 0; // Because all the edges are off
        sum_in[i] = 0;
        sum_net_out[i]=0;
        sum_net_in[i]=0;
    }
    
#ifdef MAXEDGES
    int k,kr;
    for (int i=0;i<nn;i++)
    {
        for (int j=0;j<nn;j++)
        {
            k = i + j*nn;   // position in the adjacency and weights matrix
        
            if (i!=j && on_out[i] && on_in[j]) // This stops edges to nodes with a zero in/out sum
            {
                A[k] = 1;               // switch on the edge
                active_edges[ne] = k;   // Add it to the end of the active list
                ne ++;                  // Increment number of edges
                sum_out[i] += W[k];
                sum_in[j] += W[k];
            }
        }
    }
    // Now calculate the netted positions
    double wnet;
    for (int i=0;i<nn;i++)
    {
        for (int j=i+1;j<nn;j++)
        {
            k = i + j*nn;   // position in the adjacency and weights matrix
            kr = j + i*nn;   // position of return edge
            // Calculate the difference between the weighted edges
            wnet=0;
            if (A[k])  wnet += W[k];
            if (A[kr]) wnet -= W[kr];
            if (A[k] && A[kr]) nepair++;    // count the pair
            
            // Depending which way the net flows, update the relevant sum
            if (wnet>0)
            {
                sum_net_out[i] += wnet;
                sum_net_in[j] += wnet;
            }
            else
            {
                sum_net_in[i] -= wnet;
                sum_net_out[j] -= wnet;
            }

        }
    }
    
    
    cout << "Max edges filled to ne=" << ne << " edges" << endl;
#endif
	
    energy=energynet=energyex=0;
    nnet=nex=0;
    for (int i=0;i<nn;i++)
    {
        if (target_out[i]>-1e-12)
            top[i] = pow(target_out[i],SCALE_FAC);
        else
            top[i] = INFINITY;
        
        if (target_in[i]>-1e-12)
            tip[i] = pow(target_in[i],SCALE_FAC);
        else
            tip[i] = INFINITY;
        
        if (target_net_out[i]>-1e-12)
        {
            if (target_net_out[i] < 1e-12)	// it's zero
                top_net[i] = top[i];
            else
                top_net[i] = pow(target_net_out[i],SCALE_FAC);
        }
        else
            top_net[i] = INFINITY; // This will always put the energynet to zero for this netting constraint
        
        if (target_net_in[i]>-1e-12)
        {
            if (target_net_in[i] < 1e-12)
                tip_net[i] = tip[i];
            else
                tip_net[i] = pow(target_net_in[i],SCALE_FAC);
        }
        else
            tip_net[i] = INFINITY; // This will always put the energynet to zero for this netting constraint
        
        // Count number who can monitor net
        if (target_net_out[i]>-1e-12 || target_net_in[i]>-1e-12)
            nnet ++;
        
        
        if (target_in[i]>1e-12) energy += pow2((target_in[i]-sum_in[i])/tip[i]);
        if (target_out[i]>1e-12) energy += pow2((target_out[i]-sum_out[i])/top[i]);
        // And now the netted positions
        if (target_in[i]>1e-12) energynet += pow2((target_net_in[i]-sum_net_in[i])/tip_net[i]);
        if (target_out[i]>1e-12) energynet += pow2((target_net_out[i]-sum_net_out[i])/top_net[i]);
        
        if (on_net_ex[i])
        {
            nex ++;
            
            if (fabs(target_ex_out[i])>1e-12)
                top_ex[i] = fabs(target_ex_out[i]);
            else
            {
                if (target_net_out[i]>1e-12)
                    top_ex[i] = 10*top_net[i];
                else
                    top_ex[i] = 10*tip_net[i];
            }

            energyex += pow2(target_ex_out[i]/top_ex[i]);
        }
    }

    cout << "Energy(0)=" << energy << endl;
}


double dense_hybrid_netted::total_energy()
{
    int k,l;
    
    // copy weights into activeW
    for (l=0;l<ne;l++)
    {
        k = active_edges[l];     // k is the position in the adjacency and weights matrices
        activeW[l] = W[k];
    }
    
    return(total_energy(activeW)); // overload, by default use the W matrix
}

double dense_hybrid_netted::total_energy_netted()
{
    // For netted constraints this becomes a little clunky
    // If you want to use this with CG then you may need a reverse edges array
    
    int i,j,k,kr,l;
    
    for (i=0; i<nn; i++)
    {
        cg_sum_net_in[i] = cg_sum_net_out[i]=0; // clear the netted sums
    }
    
    // Calculate all the node in/out sums
    for (l=0;l<ne;l++)
    {
        k = active_edges[l];     // k is the position in the adjacency and weights matrices
        
        // Resolve this into the from (i) and to edge (j)
        j = k/nn;
        i = k - j*nn;
        kr = j + i*nn;
        
        double wnet = W[k] - A[kr]*W[kr];
        if (wnet>0)
        {
            cg_sum_net_out[i] += wnet;
            cg_sum_net_in[j] += wnet;
        }
        
    }
    
    // Now total energy
    double tenergy=0;
    for (i=0;i<nn;i++)
    {
        if (target_net_in[i]>-1e-12 && target_in[i]>1e-12)
            tenergy += pow2((target_net_in[i] - cg_sum_net_in[i])/tip_net[i]);
        if (target_net_out[i]>-1e-12 && target_out[i]>1e-12)
            tenergy += pow2((target_net_out[i] - cg_sum_net_out[i])/top_net[i]);
    }
    return(tenergy);
}

double dense_hybrid_netted::total_energy_ex()
{
    double tenergy=0;
    int k,kr;
    for (int i=1;i<nn;i++)
    {
        if (on_net_ex[i])
        {
            k = i;
            kr = i*nn;
        
            double wnet = A[k]*W[k]-A[kr]*W[kr];
            tenergy += pow2((target_ex_out[i] - wnet)/top_ex[i]);
        }
    }
    return tenergy;
}

bool dense_hybrid_netted::rowcol_iterate()
{
    // actW is a smaller array with only active weights
    
    int i,j,k, l;
    
    double *sumhold;
    double energyi,energyj;
    
    sumhold = new double [nn];
    
    // copy weights into activeW
    for (l=0;l<ne;l++)
    {
        k = active_edges[l];     // k is the position in the adjacency and weights matrices
        activeW[l] = W[k];
    }
    
    
    int iter=0;
    bool success=false;
    double lastdeltaE, deltaE;
    int nconsec=0;
    
    while (iter < 5000 && !success)
    {
        // ROWS
        // Clear the sums
        for (int isum=0;isum<nn;isum++)
            sumhold[isum]=0;
        
        // Work out the row sum
        for (l=0;l<ne;l++)
        {
            k = active_edges[l];     // k is the position in the adjacency and weights matrices
            // Resolve this into the from (i) and to edge (j)
            j = k/nn;
            i = k - j*nn;
            
            sumhold[i] += activeW[l];
        }
        
        // Now adjust the sums
        for (l=0;l<ne;l++)
        {
            k = active_edges[l];
            j = k/nn;
            i = k - j*nn;
            
            if (target_out[i]>-1e-12)
                activeW[l] *= target_out[i]/sumhold[i];
        }
        
        energyi = total_energy(activeW);
        
        // COLUMS
        // Clear the sums
        for (int isum=0;isum<nn;isum++)
            sumhold[isum]=0;
        
        // Work out the row sum
        for (l=0;l<ne;l++)
        {
            k = active_edges[l];     // k is the position in the adjacency and weights matrices
            // Resolve this into the from (i) and to edge (j)
            j = k/nn;
            
            sumhold[j] += activeW[l];
        }
        
        // Now adjust the sums
        for (l=0;l<ne;l++)
        {
            k = active_edges[l];
            j = k/nn;
            
            if (target_in[j]>-1e-12)
                activeW[l] *= target_in[j]/sumhold[j];
        }
        
        energyj = total_energy(activeW);
        
        
        deltaE = fabs(energyj-energyi);
        
        if (deltaE < 1e-18)
        {
            success=true;
            cout << "ROW-COL-ITER: SUCCESS " << iter << endl;
            break;
        }
        
        if (iter>5)
        {
            if (fabs(1-deltaE/lastdeltaE) < 1e-14)
                nconsec++;
            else
                nconsec=0;
        }
        
        if (nconsec>5)
        {
            cout << "ROW-COL-ITER: FAIL " << iter << endl;
            success=false;
            break;
        }
        
        lastdeltaE = deltaE;
        
        iter ++; // increment counter
    }
    
    if (iter==5000)
        cout << "ROW-COL-ITER: FAIL " << iter << endl;
    
    delete sumhold;
    
    return(success);
}




double dense_hybrid_netted::total_energy(double *actW)
{
    // activeW is a smaller
    
    int i,j,k,l;
    
    for (l=0; l<nn; l++)
    {
        cg_sum_in[l] = cg_sum_out[l]=0; // clear the sums
        cg_sum_net_in[l] = cg_sum_net_out[l]=0; // clear the netted sums
    }
    
    // Calculate all the node in/out sums
    for (l=0;l<ne;l++)
    {
        k = active_edges[l];     // k is the position in the adjacency and weights matrices
        
        // Resolve this into the from (i) and to edge (j)
        j = k/nn;
        i = k - j*nn;
        
        cg_sum_out[i] += actW[l];  // actW is a smaller array with only active weights
        cg_sum_in[j] += actW[l];
    }
    
    // Now total energy
    double tenergy=0;
    for (i=0;i<nn;i++)
    {
        if (target_in[i]>1e-12)
            tenergy += pow2((target_in[i] - cg_sum_in[i])/tip[i]);
        if (target_out[i]>1e-12)
            tenergy += pow2((target_out[i] - cg_sum_out[i])/top[i]);
    }
    
    return(tenergy);
    
}

void dense_hybrid_netted::link_active_reverse()
{
    // This is a very slow function so think twice before using it a lot
    // In it's current form it scales like ne^2. That's *edges* squared.
    int i,j,k,kr, ir;
    
    for (int iact=0;iact<ne;iact++)
    {
        k = active_edges[iact];
        // Resolve this into the from (i) and to edge (j)
        j = k/nn;
        i = k - j*nn;
        kr = j + i*nn;
        
        if (A[kr])
        {
            ir=0;
            while (active_edges[ir]!=kr && ir<ne)
                ir ++;
            
            if (ir>ne-1)
            {
                cout << "link_active_reverse: Error finding reverse edge" << endl; exit(1);
            }
            
            reverse_edges[iact]=ir;
        }
        else
        {
            // No edge coming the other way
            reverse_edges[iact]=-1;
        }
        
    }
}
