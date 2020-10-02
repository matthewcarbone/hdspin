#include "Obs/base.h"
#include "Obs/psi.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

PsiConfig::PsiConfig(const FileNames fnames, const RuntimeParameters rtp) :
    Base(fnames), rtp(rtp)

{
    max_counter = (long long) log2l(ipow(10, rtp.log_N_timesteps));

    // Give the max counter a lot of space
    max_counter += 10;

    for (int ii=0; ii<max_counter; ii++){counter.push_back(0);}
    for (int ii=0; ii<max_counter; ii++){counter_IS.push_back(0);}
}


void PsiConfig::_help_step(const bool inherent_structure)
{
    long long key;

    // Use a local variable to hold the current state of the waiting time,
    // and reset it immediately since one way or another it needs to be done.
    long double local_waiting_time;

    if (inherent_structure)
    {
        local_waiting_time = waiting_time_IS;
        waiting_time_IS = 0.0;
    }
    else
    {
        local_waiting_time = waiting_time;
        waiting_time = 0.0;
    }

    // If the waiting time is <= 1, round it to 1.
    if (local_waiting_time <= 1.0){key = 0;}

    else
    {
        const long double log_t = log2l(local_waiting_time);
        key = (long long) roundl(log_t);

        // We ignore any crazy waiting times produced near the end of the
        // Gillespie dynamics since they can be chalked up to edge effects.
        if (key > max_counter - 1){return;}
    }

    // Update
    if (inherent_structure){counter_IS[key] += 1;}
    else{counter[key] += 1;}
}


void PsiConfig::step_(const long double current_waiting_time, const Vals prev,
    const Vals curr)
{
    // No matter what, we update the internal states. If the configs are logged
    // to the counters they are reset in the helper functions.
    waiting_time += current_waiting_time;
    waiting_time_IS += current_waiting_time;
    if (curr.int_rep != prev.int_rep){_help_step(false);}
    if (curr.int_rep_IS != prev.int_rep_IS){_help_step(true);}
}


// Upon destruction, automatically save to disk
PsiConfig::~PsiConfig()
{
    outfile = fopen(fnames.psi_config.c_str(), "w");
    for (int ii=0; ii<max_counter; ii++)
    {
        fprintf(outfile, "%i %lli %lli\n", ii, counter[ii], counter_IS[ii]);
    }
    // Outfile closed in base class destructor
}
