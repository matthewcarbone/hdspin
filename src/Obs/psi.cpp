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


PsiConfig::PsiConfig(const FileNames fnames, const RuntimeParameters rtp) :
    Base(fnames), rtp(rtp)
{
    max_counter = (long long) log2l(ipow(10, rtp.log_N_timesteps));

    // Give the max counter a lot of space
    max_counter += 10;

    for (int ii=0; ii<max_counter; ii++){counter.push_back(0);}
    for (int ii=0; ii<max_counter; ii++){counter_IS.push_back(0);}
}


void PsiConfig::_help_step_(const bool inherent_structure)
{

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

    const long long key = _get_key(local_waiting_time);

    // We ignore any crazy waiting times produced near the end of the
    // Gillespie dynamics since they can be chalked up to edge effects.
    if (key > max_counter - 1){return;}

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
    if (curr.int_rep != prev.int_rep){_help_step_(false);}
    if (curr.int_rep_IS != prev.int_rep_IS){_help_step_(true);}
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


PsiBasin::PsiBasin(const FileNames fnames, const RuntimeParameters rtp) :
    Base(fnames), rtp(rtp)
{
    max_counter = (long long) log2l(ipow(10, rtp.log_N_timesteps));

    // Give the max counter a lot of space
    max_counter += 10;

    for (int ii=0; ii<max_counter; ii++){counter_E.push_back(0);}
    for (int ii=0; ii<max_counter; ii++){counter_E_IS.push_back(0);}
    for (int ii=0; ii<max_counter; ii++){counter_S.push_back(0);}
    for (int ii=0; ii<max_counter; ii++){counter_S_IS.push_back(0);}
    for (int ii=0; ii<max_counter; ii++){
        counter_unique_configs_per_basin_E.push_back(0);
    }
    for (int ii=0; ii<max_counter; ii++){
        counter_unique_configs_per_basin_E_IS.push_back(0);
    }
    for (int ii=0; ii<max_counter; ii++){
        counter_unique_configs_per_basin_S.push_back(0);
    }
    for (int ii=0; ii<max_counter; ii++){
        counter_unique_configs_per_basin_S_IS.push_back(0);
    }
}


void PsiBasin::_step_E_(const long double current_waiting_time,
    const Vals prev, const Vals curr)
{

    if (prev.energy < rtp.energetic_threshold)
    {
        waiting_time_E += current_waiting_time;
        tmp_unique_configs_in_basin_E.push_back(prev.int_rep);

        if (curr.energy >= rtp.energetic_threshold)
        {
            const long double local_waiting_time = waiting_time_E;
            waiting_time_E = 0.0;

            const long long key = _get_key(local_waiting_time);
            if (key <= max_counter - 1)
            {
                counter_E[key] += 1;  // Update
            }
            
            // We also count the number of unique configs per basin
            const long double n_unique = 
                std::set<long long>(tmp_unique_configs_in_basin_E.begin(),
                    tmp_unique_configs_in_basin_E.end()).size();
            const long long key2 = _get_key(n_unique);
            if (key2 <= max_counter - 1)
            {
                counter_unique_configs_per_basin_E[key2] += 1;
            }
            tmp_unique_configs_in_basin_E.clear();
        }
    }
}


void PsiBasin::_step_E_IS_(const long double current_waiting_time,
    const Vals prev, const Vals curr)
{

    if (prev.energy_IS < rtp.energetic_threshold)
    {
        waiting_time_E_IS += current_waiting_time;
        tmp_unique_configs_in_basin_E_IS.push_back(prev.int_rep_IS);

        if (curr.energy_IS >= rtp.energetic_threshold)
        {
            const long double local_waiting_time = waiting_time_E_IS;
            waiting_time_E_IS = 0.0;

            const long long key = _get_key(local_waiting_time);
            if (key <= max_counter - 1)
            {
                counter_E_IS[key] += 1;  // Update
            }

            // We also count the number of unique configs per basin
            const long double n_unique = 
                std::set<long long>(tmp_unique_configs_in_basin_E_IS.begin(),
                    tmp_unique_configs_in_basin_E_IS.end()).size();
            const long long key2 = _get_key(n_unique);
            if (key2 <= max_counter - 1)
            {
                counter_unique_configs_per_basin_E_IS[key2] += 1;
            }
            tmp_unique_configs_in_basin_E_IS.clear();
        }
    }
}


void PsiBasin::_step_S_(const long double current_waiting_time,
    const Vals prev, const Vals curr)
{

    if (prev.energy < rtp.entropic_attractor)
    {
        waiting_time_S += current_waiting_time;
        tmp_unique_configs_in_basin_S.push_back(prev.int_rep);

        if (curr.energy >= rtp.entropic_attractor)
        {
            const long double local_waiting_time = waiting_time_S;
            waiting_time_S = 0.0;

            const long long key = _get_key(local_waiting_time);
            if (key <= max_counter - 1)
            {
                counter_S[key] += 1;  // Update
            }

            // We also count the number of unique configs per basin
            const long double n_unique = 
                std::set<long long>(tmp_unique_configs_in_basin_S.begin(),
                    tmp_unique_configs_in_basin_S.end()).size();
            const long long key2 = _get_key(n_unique);
            if (key2 <= max_counter - 1)
            {
                counter_unique_configs_per_basin_S[key2] += 1;
            }
            tmp_unique_configs_in_basin_S.clear();
        }
    }
}


void PsiBasin::_step_S_IS_(const long double current_waiting_time,
    const Vals prev, const Vals curr)
{

    if (prev.energy_IS < rtp.entropic_attractor)
    {
        waiting_time_S_IS += current_waiting_time;
        tmp_unique_configs_in_basin_S_IS.push_back(prev.int_rep_IS);

        if (curr.energy_IS >= rtp.entropic_attractor)
        {
            const long double local_waiting_time = waiting_time_S_IS;
            waiting_time_S_IS = 0.0;

            const long long key = _get_key(local_waiting_time);
            if (key <= max_counter - 1)
            {
                counter_S_IS[key] += 1;  // Update
            }

            // We also count the number of unique configs per basin
            const long double n_unique = 
                std::set<long long>(tmp_unique_configs_in_basin_S_IS.begin(),
                    tmp_unique_configs_in_basin_S_IS.end()).size();
            const long long key2 = _get_key(n_unique);
            if (key2 <= max_counter - 1)
            {
                counter_unique_configs_per_basin_S_IS[key2] += 1;
            }
            tmp_unique_configs_in_basin_S_IS.clear();
        }
    }
}


void PsiBasin::step_(const long double current_waiting_time, const Vals prev,
    const Vals curr)
{
    _step_E_(current_waiting_time, prev, curr);
    _step_S_(current_waiting_time, prev, curr);
    _step_E_IS_(current_waiting_time, prev, curr);
    _step_S_IS_(current_waiting_time, prev, curr);
}

PsiBasin::~PsiBasin()
{
    outfile = fopen(fnames.psi_basin.c_str(), "w");
    for (int ii=0; ii<max_counter; ii++)
    {
        fprintf(outfile, "%i %lli %lli %lli %lli %lli %lli %lli %lli\n",
            ii, counter_E[ii], counter_E_IS[ii], counter_S[ii],
            counter_S_IS[ii], counter_unique_configs_per_basin_E[ii],
            counter_unique_configs_per_basin_E_IS[ii],
            counter_unique_configs_per_basin_S[ii],
            counter_unique_configs_per_basin_S_IS[ii]);
    }
}
