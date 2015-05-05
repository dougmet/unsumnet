
dense_hybrid *gdh;

dense_hybrid::dense_hybrid(int nn_in, int target_ne_in)
{
    ///////////////////////////
    //// CREATE THE ARRAYS ////
    ///////////////////////////
    
    nn = nn_in;
    target_ne = target_ne_in;
    
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
    target_out = new double [nn];
    target_in = new double [nn];
    sum_out = new double [nn];
    sum_in = new double [nn];
    
    on_out = new bool [nn];
    on_in = new bool [nn];
    top = new double [nn];
    tip = new double [nn];
    
    // Conjugate gradient arrays
    activeW= new double [target_ne*2];
    gradW = new double [target_ne*2];
    cg_g = new double [target_ne*2];
    cg_h = new double [target_ne*2];
    cg_sum_out = new double [nn];
    cg_sum_in = new double [nn];
    
    // how many bloody arrays could it need?
    pcom = new double [target_ne*2];
    xicom = new double [target_ne*2];
    
    ncom = target_ne*2;
    
    gdh=this;
}


int dense_hybrid::runjob(int ncycles,
                         long  mct_schedule,
                         long  hot_time,
                         double beta0,
                         double betamax,
                         double mu,
                         double cooling_rate,
                         long max_time,
                         double cgmax, // when to attempt conjugate gradient
                         double CG_TARGET)
{
    // This function is the main job controller

    long i,j,k,l,iact;                          // general integers
    int newi,newj,nrn,newk;                     // Edge moves
    int ir, mr;                                 // for choosing moves
    long mct, max_time;                         // Monte Carlo time
    double deltaE, deltaM, w, dw, r, prob, cg_energy;
    double tenergy, temax, localmaxw;    // Working numbers

	double oldenergy;
	int nsmallderiv=0;
	
    bool targetsread=false;
	
	char fname[100];
	int nsolutions=0;
    
    
    ofstream afile, wfile, sfile, monfile;
    ifstream targetfile;

    
    sprintf(fname,"adjacency%3.3d.dat",nsolutions);
    targetfile.open(fname);
    while (targetfile.good())
    {
        nsolutions++;
        targetfile.close();
        targetfile.clear();
        sprintf(fname,"adjacency%3.3d.dat",nsolutions);
        targetfile.open(fname);
    }
    if (nsolutions>0)
        cout << "Picking up from solution " << nsolutions-1 << endl;
    


    
	// Load the input file (to be replaced with something more sophisticated)
    SCALE_FAC=1.0;
    
    beta=beta0;
    
    cout << "Checkpoint time = " << mct_schedule << endl;
    cout << "Cooling rate = " << cooling_rate << " per checkpoint time" << endl;
    cout << "Hot time = " << hot_time << " MC sweeps" << endl;
    cout << "Beta=" << beta << ", Beta max = " << betamax << endl;
    cout << endl;
    
    

    
    
    ///////////////////////////
    ///// INITIALISATION //////
    ///////////////////////////
    
    
#ifdef ENERGYBIAS
    // Collection matrix
    int nbins=10000;
    tm = new tmcontinuous(0.0, 10.0, 15, nbins, TM_NO_END_CORRECTION); // note -- you are doing end correction
#endif
		
    // Now set the targets (these will be loaded in eventually)
    
    targetfile.open("targets.in");
    if (targetfile.good())
    {
        maxw=0;
        for (i=0;i<nn;i++)
        {
            targetfile >> target_out[i];
            targetfile >> target_in[i];
            
            on_out[i] = on_in[i] = true;
            if (target_out[i]<=1e-12)
                on_out[i]=false;
            if (target_in[i]<=1e-12)
                on_in[i]=false;
            
            if (target_in[i]>maxw)
                maxw = target_in[i];
            
            if (target_out[i]>maxw)
                maxw = target_out[i];
            
            if (!targetfile.good())
                break;
        }
        
        scalew = maxw;
        maxw=1.2;
        minw=0.00001/scalew;
        
        targetsread = true;
        
        cout << "Max allowed edge weight = " << scalew*maxw << ". Largest target=" << scalew << endl;
        
        // Now scale the whole problem
        for (i=0;i<nn;i++)
        {
            target_out[i] /= scalew;
            target_in[i] /= scalew;
        }
    }
    else
    {
        cout << "Problem with targets.in" << endl; exit(1);
    }
    
	reset_arrays(); // this function zeroes A, resets W and calculates energy

    
    ////////// MOVES ////////////
	// Setup the moves we'll use
    const int Nmoves = 5;
    int moves_perMC, sum_moves;
    move_class *move[Nmoves];
    
    // Initialise the moves
    // Arguments are: Name, Moves per sweep, initial step, minimum step, maximum step, target acceptance rate
    move[0] = new move_class("Edge insertion", target_ne, 0, 0, 0, 1);       // no step size needed here
    move[1] = new move_class("Edge deletion", target_ne, 0, 0, 0, 1);        // no step size needed here
	move[2] = new move_class("Edge tweak", target_ne*4, 1.0, 0.0000001, 10.0, 0.40);   // Aiming for 40%
    move[3] = new move_class("Edge out rewire", target_ne, 0, 0, 0, 1);        // no step size needed here
    move[4] = new move_class("Edge in rewire", target_ne, 0, 0, 0, 1);        // no step size needed here
#ifdef MAXEDGES    // In the max edges run we always switch off insertions/deletions and edge moves
    move[0]->NperMC=move[1]->NperMC=move[3]->NperMC=move[4]->NperMC=0;
#endif
    
    
    // Clear some files
    monfile.open("energy.dat"); monfile.close();
    monfile.open("beta.dat"); monfile.close();
    monfile.open("step.dat"); monfile.close();
    monfile.open("success.dat"); monfile.close();
    monfile.open("ne.dat"); monfile.close();
    monfile.open("edge_corr.dat"); monfile.close();
    monfile.open("maxlocal.dat"); monfile.close();
    
    
    ////////////////////////////////////
    //////////// MAIN LOOP /////////////
    ////////////////////////////////////
    
    
    for (mct=0;mct<ncycles;mct++)
    {
#ifndef MAXEDGES    // In the max edges run we always switch off insertions/deletions and edge moves
        if (!goforquench) // restore insertions/deletions
            move[0]->NperMC=move[1]->NperMC=move[3]->NperMC=move[4]->NperMC=target_ne;
#else
        goforquench=true;
#endif
        
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
#ifdef NORETURN
                        if ((A[j + i*nn] == 0) && (A[i + j*nn] == 0))
#else
                        if (A[k] == 0)
#endif
                        {
                            // Calculate the change at the out & in nodes
                            w = W[k];

                            deltaE = w*(w - 2*(target_out[i] - sum_out[i]))/pow(top[i],2);    // Out node
                            deltaE += w*(w - 2*(target_in[j] - sum_in[j]))/pow(tip[j],2);     // In node
                            
                            // Now the change in the measure of sparsity
                            deltaM = 1 - 2*(target_ne - ne);
                            
                            // Accept or reject
                            prob = exp(-beta*(deltaE + mu*deltaM)) * vol / ((double) ne + 1.0);
#ifdef ENERGYBIAS
                            // Add to the collection matrix (also returns weight difference)
                            bias = tm->update_collection_matrix(energy, energy + deltaE, prob);
                            
							// Now apply the weight
                            if (usebias)
                                prob = exp(-beta*(deltaE + mu*deltaM) + bias ) * vol / ((double) ne + 1.0);
#endif
							
                            if (mt() < prob)
                            {
                                // Success!
                                A[k] = 1;                // switch on the edge
                                
                                active_edges[ne] = k;    // Add it to the end of the active list
                                ne ++;                          // Increment number of edges
                                
                                sum_out[i] += w; // update the totals
                                sum_in[j] += w;
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
#ifdef ENERGYBIAS
                            // If you proposed a move out of range you still need to update the collection matrix
                            tm->update_collection_matrix(energy, energy,0);
#endif
                            move[mr]->update_attempts(0); // move unsucessful
                        }
                    }
#ifdef ENERGYBIAS
                    else
                    {
                        // If you proposed a move out of range you still need to update the collection matrix
                        tm->update_collection_matrix(energy, energy,0);
                    }
#endif
                    
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
                        
#ifdef ENERGYBIAS
                        // Add to the collection matrix (also returns weight difference)
						bias = tm->update_collection_matrix(energy, energy + deltaE, prob);
						
						// Now apply the weight
                        if (usebias)
                            prob = exp(-beta*(deltaE + mu*deltaM) + bias) * ((double) ne) / vol;
#endif
                        
                        if (mt() < prob)
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
                        
                        r = move[mr]->step*(1.0 - 2*mt());    // evenly distributed on +/- stepsize
                        
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
#ifdef ENERGYBIAS
                            // Add to the collection matrix (also returns weight difference)
                            bias = tm->update_collection_matrix(energy, energy + deltaE, prob);
                            
							// Now apply the weight
                            if (usebias)
                                prob = exp(-deltaE*beta + bias);
#endif
                            
                            if (mt() < prob)
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
#ifdef ENERGYBIAS
                            // If you proposed a move out of range you still need to update the collection matrix
                            tm->update_collection_matrix(energy, energy,0);
#endif
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
#ifdef ENERGYBIAS
                            // Add to the collection matrix (also returns weight difference)
                            bias = tm->update_collection_matrix(energy, energy + deltaE, prob);
                            
							// Now apply the weight
                            if (usebias)
                                prob = exp(-deltaE*beta + bias);
#endif

                            
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
#ifdef ENERGYBIAS
                            // If you proposed a move out of range you still need to update the collection matrix
                            tm->update_collection_matrix(energy, energy,0);
#endif
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
                            
#ifdef ENERGYBIAS
                            // Add to the collection matrix (also returns weight difference)
                            bias = tm->update_collection_matrix(energy, energy + deltaE, prob);
                            
							// Now apply the weight
                            if (usebias)
                                prob = exp(-deltaE*beta + bias);
#endif
                            
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
#ifdef ENERGYBIAS
                            // If you proposed a move out of range you still need to update the collection matrix
                            tm->update_collection_matrix(energy, energy,0);
#endif
                        }
                        
                    }
                    
                }
                    break;

                    
            } // END SWITCH
            
            
        }
        //// END OF MONTE CARLO SWEEP ////

        
        
         ///////////////////////////////////////////////////////////
                        //// MEASUREMENTS ////
         ///////////////////////////////////////////////////////////
        
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
            for (mr=0;mr<Nmoves;mr++)
				move[mr]->update_step();
            
			if ((mct%(10*mct_schedule))==0)
            {
                cout <<  "  t=" << mct << ", beta=" << beta << ", step=" << move[2]->step << ":" << move[2]->success_rate << ", E=" << energy/nn << ", ne=" << ne << ":" << move[0]->success_rate << ":" << move[1]->success_rate << " rw" << move[3]->success_rate << ":" << move[4]->success_rate << endl;
                
            
                //cout << total_energy() - energy << endl;
                energy = total_energy();

                if (mct>1000)
                {
                    write_monitor("energy.dat", mct, energy/nn);
                    write_monitor("beta.dat", mct, beta);
                    write_monitor("step.dat", mct, move[2]->step);
                    write_monitor("success.dat", mct, move[2]->success_rate);
                    write_monitor("ne.dat", mct, (double) ne);
                }
				
				if ((mct%(10*mct_schedule))==0)
                {
                    afile.open("adjacency.dat");
                    for (i=0;i<nn;i++)
                    {
                        for (j=0;j<nn;j++)
                            afile << A[i + j*nn] << " ";
                        afile << endl;
                    }
                    afile.close();
                    
                    
                    wfile.open("weights.dat");
                    wfile.precision(16);
                    for (i=0;i<nn;i++)
                    {
                        for (j=0;j<nn;j++)
                            wfile << scalew*W[i + j*nn] << " ";
                        wfile << endl;
                    }
                    wfile.close();
                    
                    sfile.open("sumtargets.dat");
                    sfile.precision(16);
                    for (i=0;i<nn;i++)
                        sfile << target_in[i] << " " << sum_in[i] << " " << target_out[i] << " " << sum_out[i] << endl;
                    sfile.close();
                    
                    sfile.open("energybeta.dat");
                    sfile.precision(16);
                    temax=0;
                    for (i=0;i<nn;i++)
                    {
                        tenergy=0;
                        if (target_in[i]>1e-12)
                            tenergy += pow((target_in[i] - sum_in[i])/tip[i],2);
                        if (target_out[i]>1e-12)
                            tenergy += pow((target_out[i] - sum_out[i])/top[i],2);

                         sfile << i << " " << beta*tenergy << endl;
                        if (tenergy>temax)
                            temax=tenergy;
                    }
                    sfile.close();
                    write_monitor("maxlocal.dat", mct, temax);
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
					cout << "End of hot_time" << endl;
					energy = total_energy();
				}
                
#ifdef ENERGYBIAS
                // Kill any collection matrix data we were gathering
                tm->reset_collection();
#endif
                
            }
            else
            {
              
#ifdef ENERGYBIAS
                // Transition matrix
                if (mct%(20*mct_schedule)==0 && nrws<10)
                {
                    cout << "getting new weights " << nrws << endl;
                    tm->save_collection("collection.dat");
                    tm->getlogP();
                    // Shift weights to zero at the ends
					tm->zerologPatend();
                    tm->save_logP("logP.dat");
                    nrws++;
                }
                // End transition matrix
#endif
   
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
					
					//if (mct>max_time || (nsmallderiv>20 && beta>10000))
					if (mct>max_time || nsmallderiv>10)
					{
                        if (mct>max_time)
                            cout << endl << endl << "Gone over max time.";
                        else
                            cout << endl << endl << "Energy plateaued.";
						cout << ".... Restarting quench again" << endl << endl;
						beta=beta0;
						mct=0;
                        nsmallderiv=0;
						reset_arrays();
                        //						move[0]->NperMC=move[1]->NperMC=target_ne;
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
                            cout << endl << endl << "Rowcol fail. Restarting quench again" << endl << endl;
                            beta=beta0;
                            mct=0;
                            reset_arrays();
                        }
                        else
                        {
                            if (ne==target_ne)
                            {
                                cout << "Switching off insertions/deletions." << endl;
                                goforquench=true;
                                // Switch off insertions/deletions, we'll stick with this config
                                move[0]->NperMC=move[1]->NperMC=0;
                                
//                                cout << "Switching energy scale factor." << endl;
//                                for (int inn=0;inn<nn;inn++)
//                                {
//                                    top[inn] = pow(target_out[inn],SCALE_FAC_DROP);
//                                    tip[inn] = pow(target_in[inn],SCALE_FAC_DROP);
//                                }
//                            
//                                cout << "Energy from " << energy << flush;
//                                energy = total_energy();
//                                cout << " to " << energy << endl;
//                                beta /= 10;
                            }
                        }
                    }

                    
                    // Kill edge rewiring when that stops
                    if (move[3]->success_rate < 1e-3 && move[3]->NperMC > 0.1)
                    {
                        if (rowcol_iterate())
                        {
                            cout << "Freezing edge swaps." << endl;
                            move[3]->NperMC=move[4]->NperMC=0;
                        }
                    }
                    
                }
                
                
                ///////////////////////////////////////////////////////////
                    /////////////////  FINAL ANSWER ///////////////
                ///////////////////////////////////////////////////////////
                
                // Conjugate gradient
                if (mct%(50*mct_schedule)==0)
                if ((energy<cgmax*nn && rowcol_iterate()) || beta>betamax)
                {
                    
                    bool hittarget=false;
                    
                    if (energy<CG_TARGET*nn)
                    {
                        hittarget=true;
                    }
                    else
                    {
                        conjugate_gradient();
                        cg_energy = total_energy(activeW);
                        cout << "cg_energy/nn=" << cg_energy/nn << endl;
                    
                        if (cg_energy<CG_TARGET*nn)
                        {
                            // check weights
                            bool allwpos=true;
                            for (l=0;l<ne;l++)
                            {
                                if (activeW[l]<minw)
                                    cout << "Small edge " << activeW[l] << endl;
                                
                                if (activeW[l]<minw*0)
                                {
                                    cout << activeW[l] << " ";
                                    allwpos=false;
                                    break;
                                }
                            }
                            
                            if (allwpos)
                            {
                                hittarget=true;
                                
                                // copy weights
                                for (l=0;l<ne;l++)
                                {
                                    k = active_edges[l];
                                    W[k] = activeW[l];
                                }
                            }
                            else
                            {
                                cout << endl << "Invalid solution. Keep cooling." << endl;
                            }
                                
                        }
                    }
                    
                    
                    if (hittarget)
                    {
                        sprintf(fname,"adjacency%3.3d.dat",nsolutions);
                        afile.open(fname);
                        for (i=0;i<nn;i++)
                        {
                            for (j=0;j<nn;j++)
                                afile << A[i + j*nn] << " ";
                            afile << endl;
                        }
                        afile.close();
                        
                        sprintf(fname,"weights%3.3d.dat",nsolutions);
                        wfile.open(fname);
                        wfile.precision(16);
                        for (i=0;i<nn;i++)
                        {
                            for (j=0;j<nn;j++)
                                wfile << scalew*W[i + j*nn] << " ";
                            wfile << endl;
                        }
                        wfile.close();
                        
                        sfile.open("sumtargets.dat");
                        sfile.precision(16);
                        for (i=0;i<nn;i++)
                            sfile << target_in[i] << " " << cg_sum_in[i] << " " << target_out[i] << " " << cg_sum_out[i] << endl;
                        sfile.close();
                        
                        cout << endl << ">>>>>> hit target " << nsolutions << " >>>>>>" << endl << endl;

                        nsolutions++;
                        beta=beta0;
                        mct=0;
                        reset_arrays();
//						move[0]->NperMC=move[1]->NperMC=target_ne;
                        
                        if (nsolutions>100)
                            exit(0);
                        
                    }

                }

				
				
            }
            

            
            
                        
            
//            if (energy<1e-5)
//            {
//                cout << "hit energy target by MC, E/n=" << energy/nn << " " << total_energy()/nn << endl;
//                exit(0);
//            }

        }
        
    } // end mct loop
    
    
	return 0;
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
    
#ifdef MAXEDGES
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
    cout << "Max edges filled to ne=" << ne << " edges" << endl;
#endif
	
    energy=0;
    for (int i=0;i<nn;i++)
    {
        top[i] = pow(target_out[i],SCALE_FAC);
        tip[i] = pow(target_in[i],SCALE_FAC);
        
        if (target_in[i]>1e-12) energy += pow((target_in[i]-sum_in[i])/tip[i],2); // initialise energy, remember all edges are off
        if (target_out[i]>1e-12) energy += pow((target_out[i]-sum_out[i])/top[i],2); // initialise energy, remember all edges are off
    }
    cout << "Energy(0)=" << energy << endl;
}

void gradenergy(double *p, double *xi)
{
    gdh->gradient_energy(p,xi);
}

double energyfunc(double *p)
{
    return(gdh->total_energy(p));
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

void dense_hybrid::gradient_energy(double *actW, double *xi)
{
    // actW is a smaller array with only active weights
    
    int i,j,k, l;
    
    // Calculate all the node in/out sums
    // DO WE NEED TO DO THIS IF WE'RE DOING IT IN THE ENERGY? CHECK minimise.cpp
    for (l=0; l<nn; l++)
        cg_sum_in[l] = cg_sum_out[l]=0; // clear the sums
    
    for (l=0;l<ne;l++)
    {
        k = active_edges[l];     // k is the position in the adjacency and weights matrices
        
        // Resolve this into the from (i) and to edge (j)
        j = k/nn;
        i = k - j*nn;
        
        cg_sum_out[i] += actW[l];
        cg_sum_in[j] += actW[l];
    }
    
    // Now work out the gradients
    for (l=0;l<ne;l++)
    {
        k = active_edges[l];     // k is the position in the adjacency and weights matrices
        
        // Resolve this into the from (i) and to edge (j)
        j = k/nn;
        i = k - j*nn;
        
        xi[l]=0;
        if (target_out[i]>1e-12)
            xi[l] = -2*((target_out[i] - cg_sum_out[i])/pow(top[i],2))*actW[l];
        if (target_in[j]>1e-12)
            xi[l] += -2*((target_in[j] - cg_sum_in[j])/pow(tip[j],2))*actW[l];
    }
}

void dense_hybrid::conjugate_gradient()
{
    long k,l,iter;
	double fret;
    
    // copy weights into activeW
    for (l=0;l<ne;l++)
    {
        k = active_edges[l];     // k is the position in the adjacency and weights matrices
        activeW[l] = W[k];
    }
    
    // Now call the conjugate gradient minimisation routines
	
	frprmn(activeW, ne, FTOL, &iter, &fret, energyfunc, gradenergy, cg_g, cg_h, gradW);
	
	cout << "Relaxed to precision FTOL=" << FTOL << " in " << iter << " iterations." << endl;
}
