/*
 
 This file is part of dense_hybrid.
 
 dense_hybrid is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 2 of the License, or
 (at your option) any later version.
 
 dense_hybrid is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with dense_hybrid.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#include "dense_hybrid.h"

dense_hybrid::dense_hybrid(int nn_in, int target_ne_in,
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
    
    // They're 1D to keep pointers simple
    // I've put them in COLUMN-MAJOR format as this is compatible with blas and lapack
    // libraries incase we want eigenvectors etc.
    // To get element A_ij address A[i + j*nn]
    A = new int [nn*nn];
    W = new double [nn*nn];
    active_edges = new long [max_ne];     // first ne elements give active edges
    active_edges0 = new long [max_ne];    // For correlation functions
    vol = nn*(nn-1);                           // total number of possible edges
    
    // These arrays store node properties
    target_out = new double [nn];
    target_in = new double [nn];
    sum_out = new double [nn];
    sum_in = new double [nn];
    
    on_out = new bool [nn];
    on_in = new bool [nn];
    top = new double [nn];
    tip = new double [nn];
    
    // Conjugate gradient arrays
    activeW= new double [max_ne];
    gradW = new double [max_ne];
    cg_g = new double [max_ne];
    cg_h = new double [max_ne];
    cg_sum_out = new double [nn];
    cg_sum_in = new double [nn];
    
    
}

dense_hybrid::~dense_hybrid()
{
    delete[] A;
    delete[] W;
    delete[] active_edges;
    delete[] active_edges0;
    delete[] target_out;
    delete[] target_in;
    delete[] sum_out;
    delete[] sum_in;
    delete[] on_out;
    delete[] on_in;
    delete[] top;
    delete[] tip;
    delete[] activeW;
    delete[] gradW;
    delete[] cg_g;
    delete[] cg_h;
    delete[] cg_sum_out;
    delete[] cg_sum_in;
}

// Using std::vectors for compatibility outside. I appreciate the jarring
// change of style. We'll copy into our old school arrays and move along
void dense_hybrid::init_targets(const std::vector<double> & target_out_i,
                                const std::vector<double> & target_in_i)
{
    // Find maximum and find which nodes have no out/in going edges
    maxw=0;
    for (int i=0;i<nn;i++)
    {
        on_out[i] = on_in[i] = true;
        if (target_out_i[i]<=1e-12)
            on_out[i]=false;
        if (target_in_i[i]<=1e-12)
            on_in[i]=false;
        
        if (target_in_i[i]>maxw)
            maxw = target_in_i[i];
        
        if (target_out_i[i]>maxw)
            maxw = target_out_i[i];
        
    }

    scalew = maxw;
    maxw=1.2;
    minw=0.00001/scalew;

    // Now scale the whole problem
    for (int i=0;i<nn;i++)
    {
        target_out[i] = target_out_i[i] / scalew;
        target_in[i] = target_in_i[i] / scalew;
    }

}

int dense_hybrid::runjob(long  mct_schedule,
                         long  hot_time,
                         double beta0,
                         double betamax,
                         double mu0,
                         double cooling_rate,
                         long max_time,
                         double cgmax)
{
    // This function is the main job controller

    long mct;                           // Monte Carlo time
    int nsmallderiv=0;
    double oldenergy=1;


    ///////////////////////////
    ///// INITIALISATION //////
    ///////////////////////////
    
    SCALE_FAC=1.0;
    
    beta=beta0;
    mu=mu0;
	reset_arrays(); // this function zeroes A, resets W and calculates energy

    
    ////////// MOVES ////////////
	// Setup the moves we'll use
    const int Nmoves = 5;
    move_class *move[Nmoves];
    
    // Initialise the moves
    // Arguments are: Name, Moves per sweep, initial step, minimum step, maximum step, target acceptance rate
    move[0] = new move_class("Edge insertion", target_ne, 0, 0, 0, 1);       // no step size needed here
    move[1] = new move_class("Edge deletion", target_ne, 0, 0, 0, 1);        // no step size needed here
	move[2] = new move_class("Edge tweak", target_ne*4, 1.0, 0.0000001, 10.0, 0.40);   // Aiming for 40%
    move[3] = new move_class("Edge out rewire", target_ne, 0, 0, 0, 1);        // no step size needed here
    move[4] = new move_class("Edge in rewire", target_ne, 0, 0, 0, 1);        // no step size needed here
if (MAXEDGES)    // In the max edges run we always switch off insertions/deletions and edge moves
    move[0]->NperMC=move[1]->NperMC=move[3]->NperMC=move[4]->NperMC=0;
    
    
    ////////////////////////////////////
    //////////// MAIN LOOP /////////////
    ////////////////////////////////////
    
    
    for (mct=0;mct<max_time;mct++)
    {
        if(!MAXEDGES)    // In the max edges run we always switch off insertions/deletions and edge moves
        {
            if (!goforquench) // restore insertions/deletions
                move[0]->NperMC=move[1]->NperMC=move[3]->NperMC=move[4]->NperMC=target_ne;
        }
        else
        {
            goforquench=true;
        }

        ////////////////////////////
        ////    MAIN MC LOOP    ////
        ////////////////////////////
        
        mc_sweep(move, Nmoves);
    
        
        ////////////////////////////
        ////    MEASUREMENTS    ////
        ////////////////////////////
        
        if ((mct%mct_schedule)==0)
        {
            // Change the temperature in the annealing scheme
            if (mct>hot_time)
            {
                if (beta<betamax)
                {
                        beta *= cooling_rate;
                }
            }
            
            // Update the step sizes based on the acceptance rate
            for (int mr=0;mr<Nmoves;mr++)
				move[mr]->update_step();
            
			if ((mct%(10*mct_schedule))==0)
            {
#ifdef BUILD_FOR_R
                if(verbose) {
                    Rcpp::Rcout <<  "  t=" << mct << ", beta=" << beta
                    << ", step=" << move[2]->step << ":" << move[2]->success_rate
                    << ", E=" << energy/nn << ", ne=" << ne << ":" << move[0]->success_rate
                    << ":" << move[1]->success_rate << " rw" << move[3]->success_rate
                    << ":" << move[4]->success_rate << endl;
                }
#endif
                
            
                //cout << total_energy() - energy << endl;
                energy = total_energy();

				
            }

            ////////////////////////////
            ////   QUENCH CONTROL   ////
            ////////////////////////////
            
            if (mct<=hot_time)
            {
                // Save a list of who is switched on before we start cooling
                ne0=ne;
                for (int i=0;i<ne;i++)
                    active_edges0[i]=active_edges[i];
				
				if (mct==hot_time)
				{
#ifdef BUILD_FOR_R
					if(verbose)
                        Rcpp::Rcout << "End of hot_time" << endl;
#endif
					energy = total_energy();
				}                
            }
            else
            {
              
                int i,k;
                for (k=0, i=0;i<ne0;i++)
                    k += A[active_edges0[i]]; // sum all edges who are still on
                //write_monitor("edge_corr.dat", mct, ((double) k)/ne0);

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
                        for(int i=0;i<Nmoves;i++) delete move[i];
                        return(DH_FAIL_PLATEAU);
					}
				}
				
                // Row-Column iteration solution finder
                if (mct%(50*mct_schedule)==0)
                {

                   // if (energy<nn)
                    if (move[0]->success_rate < 1e-5 && !goforquench)
                    {
                        if (!rowcol_iterate())
                        {
                            // This structure cannot work
                            return(DH_FAIL_ROWCOL);
                        }
                        else
                        {
                            if (ne==target_ne)
                            {
#ifdef BUILD_FOR_R
                                if(verbose)
                                    Rcpp::Rcout << "Switching off insertions/deletions." << endl;
#endif
                                goforquench=true;
                                // Switch off insertions/deletions, we'll stick with this config
                                move[0]->NperMC=move[1]->NperMC=0;
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
                                Rcpp::Rcout << "Freezing edge swaps." << endl;
#endif
                            move[3]->NperMC=move[4]->NperMC=0;
                        }
                    }
                    
                }
                
                
                ////////////////////////////
                ////    FINAL ANSWER?   ////
                ////////////////////////////
                
                if (mct%(50*mct_schedule)==0 && energy<cgmax*nn)
                {
                    if (rowcol_iterate())
                    {
                        for(int i=0;i<Nmoves;i++) delete move[i];
                        return(DH_SUCCESS);
                    }
                }

            } // hot time if

#ifdef BUILD_FOR_R
            // Check in with R
            Rcpp::checkUserInterrupt();
#endif
        } // mct_schedule if
        
    } // end mct loop
    
    // Hit the buffers without success
    for(int i=0;i<Nmoves;i++) delete move[i];
	return (DH_FAIL_TIME_OUT);
}


//////////////////////////////////
//      Monte Carlo Sweep       //
//////////////////////////////////

int dense_hybrid::mc_sweep(move_class ** move, int Nmoves)
{
    
    int moves_perMC, sum_moves;
    long i,j,k,iact;                  // general integers
    int newi,newj,nrn,newk;             // Edge moves
    int ir, mr;                         // for choosing moves
    double deltaE, deltaM, w, dw;
    double r, prob;
    double localmaxw;   // Working numbers
    
    
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
                
                
                // We need to check if the edge and the edge going the other way is on
                // Not allowing loops length 2 at the moment.
                //if (A[k] == 0)
                
                if (on_out[i] && on_in[j]) // This stops edges to nodes with a zero in/out sum
                {
                    int move_succeed=0;
                    
                    if (A[k] == 0 && ne < max_ne)
                    {
                        if (A[j + i*nn] == 0 || !NORETURN)
                        {
                            // Calculate the change at the out & in nodes
                            w = W[k];
                            
                            deltaE = w*(w - 2*(target_out[i] - sum_out[i]))/pow(top[i],2);    // Out node
                            deltaE += w*(w - 2*(target_in[j] - sum_in[j]))/pow(tip[j],2);     // In node
                            
                            // Now the change in the measure of sparsity
                            deltaM = 1 - 2*(target_ne - ne);
                            
                            // Accept or reject
                            prob = exp(-beta*(deltaE + mu*deltaM)) * vol / ((double) ne + 1.0);
                            
                            if (mt.rand53() < prob)
                            {
                                // Success!
                                A[k] = 1;                // switch on the edge
                                
                                active_edges[ne] = k;    // Add it to the end of the active list
                                ne ++;                          // Increment number of edges
                                
                                sum_out[i] += w; // update the totals
                                sum_in[j] += w;
                                energy += deltaE;
                                
                                move_succeed=1;
                                
                                
                            }
                        }
                    }
                    
                    move[mr]->update_attempts(move_succeed); // move sucessful or not
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
                    
                    // What's the edge weight for the chosen edge
                    w = W[k];
                    
                    // Calculate how this will change the energies
                    deltaE = w*(w + 2*(target_out[i] - sum_out[i]))/pow(top[i],2);    // Out node
                    deltaE += w*(w + 2*(target_in[j] - sum_in[j]))/pow(tip[j],2);     // In node
                    
                    // Now the change in the measure of sparsity
                    deltaM = 1 + 2*(target_ne - ne);
                    
                    
                    // Accept or reject
                    prob = exp(-beta*(deltaE + mu*deltaM)) * ((double) ne) / vol;
                    
                    
                    if (mt.rand53() < prob)
                    {
                        
                        // Success!
                        A[k] = 0;                               // Switch off the edge
                        active_edges[iact]=active_edges[ne-1];  // Replace the entry in the active list with the one at the end
                        ne --;                                  // Decrement the number of edges
                        
                        sum_out[i] -= w; // update the totals
                        sum_in[j] -= w;
                        
                        energy += deltaE;
                        
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
                    
                    r = move[mr]->step*(1.0 - 2*mt.rand53());    // evenly distributed on +/- stepsize
                    
                    if (target_out[i]>target_in[j])
                        localmaxw=1.2*target_in[j];
                    else
                        localmaxw=1.2*target_out[i];
                    
                    dw = r*localmaxw;
                    
                    //                        if (mt.randInt(1))
                    //                            dw = r*target_in[j];                  // in W
                    //                        else
                    //                            dw = r*target_out[i];                  // in W
                    
                    
                    //r = exp(r);                         // random walk in log(W)
                    //dw = W[k]*(r-1);                    // random walk in W
                    
                    
                    
                    //if ((W[k] + dw > minw) && (W[k] + dw < maxw))
                    if ((W[k] + dw > minw) && (W[k] + dw < localmaxw))
                    {
                        deltaE = dw*(dw - 2*(target_out[i] - sum_out[i]))/pow(top[i],2);    // Out node
                        deltaE += dw*(dw - 2*(target_in[j] - sum_in[j]))/pow(tip[j],2);     // In node
                        
                        // Accept or reject
                        prob = exp(-deltaE*beta);
                        
                        
                        if (mt.rand53() < prob)
                        {
                            W[k] += dw;
                            sum_out[i] += dw;
                            sum_in[j] += dw;
                            energy += deltaE;
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
                    
                    if (nrn<nn && w<localmaxw)
                    {
                        
                        deltaE = w*(w - 2*(target_in[newj] - sum_in[newj]))/pow(tip[newj],2);   // New in node
                        deltaE += w*(w + 2*(target_in[j] - sum_in[j]))/pow(tip[j],2);           // Old in node
                        
                        // Accept or reject
                        prob = exp(-deltaE*beta);
                        
                        if (mt.rand53() < prob)
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
                    
                    if (nrn<nn && w<localmaxw)
                    {
                        deltaE = w*(w + 2*(target_out[i] - sum_out[i]))/pow(top[i],2);          // Old out node
                        deltaE += w*(w - 2*(target_out[newi] - sum_out[newi]))/pow(top[newi],2);// New out node
                        
                        // Accept or reject
                        prob = exp(-deltaE*beta);
                        
                        
                        if (mt.rand53() < prob)
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
                
                
        } // END SWITCH
        
        
    }
    //// END OF MONTE CARLO SWEEP ////
    
    return(0); // success
    
}


//////////////////////////////////
// Conjugate gradient functions
//////////////////////////////////

void dense_hybrid::reset_arrays()
{
    goforquench=false;
    
	// Initialisation of edges
    for (int i=0;i<nn*nn;i++)
    {
        A[i]=0;         // All switched OFF
        W[i]=minw;   // All weights equal (could set to anything)
    }
    ne=0;               // No edges at the start
    
    // Initialise the sums of edges
    for (int i=0;i<nn;i++)
    {
        sum_out[i] = 0; // Because all the edges are off
        sum_in[i] = 0;
    }
    
    if (MAXEDGES)
    {
        for (int i=0;i<nn;i++)
        {
            for (int j=0;j<nn;j++)
            {
                int k = i + j*nn;   // position in the adjacency and weights matrix
            
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
//        cout << "Max edges filled to ne=" << ne << " edges" << endl;
    }
	
    energy=0;
    for (int i=0;i<nn;i++)
    {
        top[i] = pow(target_out[i],SCALE_FAC);
        tip[i] = pow(target_in[i],SCALE_FAC);
        
        if (target_in[i]>1e-12) energy += pow((target_in[i]-sum_in[i])/tip[i],2); // initialise energy, remember all edges are off
        if (target_out[i]>1e-12) energy += pow((target_out[i]-sum_out[i])/top[i],2); // initialise energy, remember all edges are off
    }
#ifdef BUILD_FOR_R
    if(verbose)
        Rcpp::Rcout << "Energy(0)=" << energy << endl;
#endif
}


double dense_hybrid::total_energy()
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

bool dense_hybrid::rowcol_iterate()
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
            
            activeW[l] *= target_out[i]/sumhold[i];
        }
        
        energyi = total_energy(activeW);
        
        // COLUMS
        // Clear the sums
        for (int isum=0;isum<nn;isum++)
            sumhold[isum]=0;
        
        // Work out the column sum
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
            
            activeW[l] *= target_in[j]/sumhold[j];
        }
        
        energyj = total_energy(activeW);
        
        
        deltaE = fabs(energyj-energyi);
        
        if (deltaE < 1e-18)
        {
            success=true;
//            cout << "ROW-COL-ITER: SUCCESS " << iter << endl;
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
//            cout << "ROW-COL-ITER: FAIL " << iter << endl;
            success=false;
            break;
        }
        
        lastdeltaE = deltaE;
        
        iter ++; // increment counter
    }
    
//    if (iter==5000)
//        cout << "ROW-COL-ITER: FAIL " << iter << endl;
    
    delete[] sumhold;
    
    return(success);
}




double dense_hybrid::total_energy(double *actW)
{
    // activeW is a smaller
    
    int i,j,k, l;
        
    for (l=0; l<nn; l++)
        cg_sum_in[l] = cg_sum_out[l]=0; // clear the sums
    
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
            tenergy += pow((target_in[i] - cg_sum_in[i])/tip[i],2);
        if (target_out[i]>1e-12)
            tenergy += pow((target_out[i] - cg_sum_out[i])/top[i],2);
    }
    
    return(tenergy);
    
}
