/*
 *  move_class.cpp
 *  
 *
 *  Created by Douglas Ashton on 08/04/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

move_class::move_class(const char *setname, int setnpermc, double setstep, double setmin, double setmax, double settarget)
{
    strcpy(name, setname);
    NperMC  = setnpermc;
    step    = setstep;
    minstep = setmin;
    maxstep = setmax;
    target  = settarget;
    Nattempt=0;
    Nsuccess=0;
    
    pid_sum=0;
    npid=0;
    
    Nops=Nopf=0;
    opsuccess=opfail=0.0;
    
    fup = 1.045;
    fdown = 1/1.05;
    
}

void move_class::update_attempts(bool outcome, int op)
{
    update_attempts(outcome, (double) op);
}


void move_class::update_attempts(bool outcome, double op)
{
    if (outcome)
    {
        opsuccess += op;
        Nops++;
    }
    else
    {
        opfail += op;
        Nopf++;
    }
    
    update_attempts(outcome); // onward
}

void move_class::update_attempts(bool outcome)
{
    Nsuccess += outcome;
    Nattempt ++;
}

void move_class::reset_attempts()
{
    Nsuccess = Nattempt = 0;
    opsuccess = opfail = 0.0;
}

void move_class::update_step()
{
    if (Nops>0)
        ops=opsuccess/Nops;
    else
        ops = 0.0;
    
    if (Nopf>0)
        opf = opfail/Nopf;
    else
        opf=0.0;
    Nops=Nopf=0;
    opsuccess=opfail=0.0;
    
    
    if (Nattempt>1000)
    {
        success_rate = ((double) Nsuccess)/Nattempt;
        if ((success_rate > target*fup) && step*fup < maxstep)
        {
            step *= fup;
         //   cout << name << ": rate " << ((double) Nsuccess)/Nattempt << ", step up to  =" << step << endl;
            Nsuccess=Nattempt=0;
        }
        else if ((success_rate < target*fdown) && step*fdown > minstep)

        {
            step *= fdown;
         //   cout << name << ": rate " << ((double) Nsuccess)/Nattempt << ", step down to=" << step << endl;
            Nsuccess=Nattempt=0;
        }
        else if (Nattempt > 1000000)
        {
            // Stop it accumulating too much data so it can't move
            Nattempt /= 100;
            Nsuccess /= 100;
        }
    }
    

}

void move_class::pid_update_step()
{
    // Update the move based on a PID controller
    double measured, error;
    
    measured = ((double) Nsuccess)/Nattempt;
    error = target - measured;  // How far off are we?
    pid_sum += error;           // Update the integral
    
    // Calculate the new step, not using derivative here
    step = pid_p * (error + pid_i * pid_sum);
    
    
    
    // Make sure we stay in bounds
    if (step > maxstep)
        step = maxstep;
    else if (step < minstep)
        step = minstep;
    

    
    npid ++;
}

