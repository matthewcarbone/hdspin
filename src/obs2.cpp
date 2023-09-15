/**
 * Two-point correlators such as psi and Pi.
 */

#include <math.h>

#include "obs2.h"

long long _get_key(const long double local_waiting_time)
{
    long long key;

    // If the waiting time is <= 1, round it to 1.
    if (local_waiting_time <= 1.0){key = 0;}

    else
    {
        const long double log_t = log2l(local_waiting_time);
        key = (long long) roundl(log_t);
    }

    return key;
}


PsiConfig::PsiConfig(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{

    spin_system_ptr = &spin_system;

    _max_counter = (long long) log2l(ipow(10, params.log10_N_timesteps));

    // Give the max counter a lot of space
    _max_counter += 10;

    for (int ii=0; ii<_max_counter; ii++){_counter.push_back(0);}
}

void PsiConfig::step(const double current_waiting_time)
{

    const parameters::StateProperties prev = spin_system_ptr->get_previous_state();
    const parameters::StateProperties curr = spin_system_ptr->get_current_state();

    // No matter what, we update the internal states. If the configs are logged
    // to the counters they are reset in the helper functions.
    _waiting_time += current_waiting_time;

    // Now we attempt to log things. If the previous configuration representation
    // is different than the current one, we simply return.
    if (curr.state == prev.state){return;}

    // Otherwise we get the key
    const long long key = _get_key(_waiting_time);

    // Reset the waiting time
    _waiting_time = 0.0;

    // We ignore any crazy waiting times produced near the end of the
    // Gillespie dynamics since they can be chalked up to edge effects.
    if (key > _max_counter - 1)
    {
        // We also keep track of the number of times this happens
        // Hopefully it's very rare
        _out_of_counter += 1;
        return;
    }

    // Update the counter with the key
    _counter[key] += 1;
}


PsiConfig::~PsiConfig()
{
    outfile = fopen(fnames.psi_config.c_str(), "w");
    for (int ii=0; ii<_max_counter; ii++)
    {
        fprintf(outfile, "%lli\n", _counter[ii]);
    }
    fclose(outfile);
}
