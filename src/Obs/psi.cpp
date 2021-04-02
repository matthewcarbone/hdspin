#include <set>

#include "Obs/base.h"
#include "Obs/psi.h"
#include "Utils/utils.h"
#include "Utils/structures.h"


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


PsiConfig::PsiConfig(const FileNames fnames, const RuntimeParameters rtp,
    const bool inherent_structure) : Base(fnames),
    inherent_structure(inherent_structure), rtp(rtp)
{
    max_counter = (long long) log2l(ipow(10, rtp.log_N_timesteps));

    // Give the max counter a lot of space
    max_counter += 10;

    for (int ii=0; ii<max_counter; ii++){counter.push_back(0);}
}


void PsiConfig::_help_step_()
{
    const long long key = _get_key(waiting_time);
    waiting_time = 0.0;

    // We ignore any crazy waiting times produced near the end of the
    // Gillespie dynamics since they can be chalked up to edge effects.
    if (key > max_counter - 1){return;}

    // Update
    counter[key] += 1;
}


void PsiConfig::step_(const long double current_waiting_time, const Vals prev,
    const Vals curr)
{
    // No matter what, we update the internal states. If the configs are logged
    // to the counters they are reset in the helper functions.
    waiting_time += current_waiting_time;
    if (inherent_structure){if (curr.int_rep != prev.int_rep){_help_step_();}}
    else{if (curr.int_rep_IS != prev.int_rep_IS){_help_step_();}}
}


// Upon destruction, automatically save to disk
PsiConfig::~PsiConfig()
{
    if (inherent_structure){outfile = fopen(fnames.psi_config.c_str(), "w");}
    else{outfile = fopen(fnames.psi_config_IS.c_str(), "w");}
    for (int ii=0; ii<max_counter; ii++)
    {
        fprintf(outfile, "%i %lli\n", ii, counter[ii]);
    }
    // Outfile closed in base class destructor
}


PsiBasin::PsiBasin(const FileNames fnames, const RuntimeParameters rtp,
    const bool inherent_structure, const bool energetic_barrier) :
    Base(fnames), inherent_structure(inherent_structure),
    energetic_barrier(energetic_barrier), rtp(rtp)
{
    max_counter = (long long) log2l(ipow(10, rtp.log_N_timesteps));

    // Give the max counter a lot of space
    max_counter += 10;

    for (int ii=0; ii<max_counter; ii++){counter.push_back(0);}
    for (int ii=0; ii<max_counter; ii++)
    {
        counter_unique_configs_per_basin.push_back(0);
    }

    // Set the threshold
    if (energetic_barrier){threshold = rtp.energetic_threshold;}
    else{threshold = rtp.entropic_attractor;}
}

void PsiBasin::step_(const long double current_waiting_time,
    const Vals prev, const Vals curr)
{

    double _prev_energy, _curr_energy;
    long long _prev_int_rep;
    if (inherent_structure)
    {
        _prev_energy = prev.energy_IS;
        _curr_energy = curr.energy_IS;
        _prev_int_rep = prev.int_rep_IS;
    }    
    else
    {
        _prev_energy = prev.energy;
        _curr_energy = curr.energy;
        _prev_int_rep = prev.int_rep;
    }

    if (_prev_energy < threshold)
    {
        waiting_time += current_waiting_time;
        tmp_unique_configs_in_basin.push_back(_prev_int_rep);

        if (_curr_energy >= threshold)
        {
            const long double local_waiting_time = waiting_time;
            waiting_time = 0.0;

            const long long key = _get_key(local_waiting_time);
            if (key <= max_counter - 1){counter[key] += 1;}

            // We also count the number of unique configs per basin
            const long double n_unique = 
                std::set<long long>(tmp_unique_configs_in_basin.begin(),
                    tmp_unique_configs_in_basin.end()).size();
            const long long key2 = _get_key(n_unique);
            if (key2 <= max_counter - 1)
            {
                counter_unique_configs_per_basin[key2] += 1;
            }
            tmp_unique_configs_in_basin.clear();
        }
    }
}


PsiBasin::~PsiBasin()
{
    if (energetic_barrier && !inherent_structure)
    {
        outfile = fopen(fnames.psi_basin_E.c_str(), "w");    
    }
    else if (!energetic_barrier && !inherent_structure)
    {
        outfile = fopen(fnames.psi_basin_S.c_str(), "w");    
    }
    else if (energetic_barrier && inherent_structure)
    {
        outfile = fopen(fnames.psi_basin_E_IS.c_str(), "w");    
    }
    else if (!energetic_barrier && inherent_structure)
    {
        outfile = fopen(fnames.psi_basin_S_IS.c_str(), "w");    
    }
    else
    {
        std::cout << "WARNING: data may not have been saved in psi basin!"
        << std::endl;
    }
    for (int ii=0; ii<max_counter; ii++)
    {
        fprintf(outfile, "%i %lli %lli\n",
            ii, counter[ii], counter_unique_configs_per_basin[ii]);
    }
}
