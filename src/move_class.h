/*
 *  move_class.h
 *  
 *
 *  Created by Douglas Ashton on 08/04/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DOUG_MOVE_CLASS_H
#define DOUG_MOVE_CLASS_H

class move_class
	{
	// The particle class can have as much or as little as you need.
	
	public:
        // Variables
        char name[50];      // Unique name
        int Nattempt;       // Attempted moves
        int Nsuccess;       // Successful moves
        int NperMC;         // How many moves per MC sweep?
        double success_rate;// Save for diagnostics
        
        double step;        // Step size
        double minstep, maxstep; // Allowed range for step size
        double target;
        
        double fup, fdown;  // Factor to raise or lower step
        
        double pid_p, pid_i, pid_d; // PID controller parameters
        double pid_sum;
        int npid;
        
        int Nops,Nopf;
        double opsuccess, opfail, ops, opf;   // Some order parameter you may want to average
        
        // Methods
        move_class(const char *setname, int setnpermc, double setstep, double setmin, double setmax, double settarget);

        void update_step();
        void pid_update_step(); // Update method based on PID controller
        void update_attempts(bool outcome);
        void update_attempts(bool outcome, double op);
        void update_attempts(bool outcome, int op);
        void reset_attempts();
        
		
	};

#include "move_class.cpp"

#endif // header protection

