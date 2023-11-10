/**
 * Two-point correlators such as psi and Pi.
 */

#include <math.h>

#include "obs2.h"


long long _get_max_counter(const parameters::SimulationParameters params)
{
    long long mc = (long long) log2l(ipow(10, params.log10_N_timesteps));
    mc += 10;  // Give the max counter a lot of space
    return mc;
}


/**
 * @brief Fills a long long vector with zeros up to _max_counter
 */
void _fill_counter(const long long _max_counter, std::vector<long long> &counter)
{
    for (int ii=0; ii<_max_counter; ii++){counter.push_back(0);}
}


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


// PSI CONFIG -----------------------------------------------------------------

PsiConfig::PsiConfig(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{
    spin_system_ptr = &spin_system;
    _max_counter = _get_max_counter(params);
    _fill_counter(_max_counter, _counter);
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


// PSI BASIN ------------------------------------------------------------------


void PsiBasin::_init_E_objects()
{
    _fill_counter(_max_counter, E_objects._counter);
    _fill_counter(_max_counter, E_objects._counter_unique_configs_per_basin);
    E_objects.threshold = params.energetic_threshold;
    E_objects.outfile = fopen(fnames.psi_basin_E.c_str(), "w");
}

void PsiBasin::_init_S_objects()
{
    _fill_counter(_max_counter, S_objects._counter);
    _fill_counter(_max_counter, S_objects._counter_unique_configs_per_basin);
    S_objects.threshold = params.entropic_attractor;
    if (S_objects.threshold_valid)
    {
        S_objects.outfile = fopen(fnames.psi_basin_S.c_str(), "w");
    }
}


PsiBasin::PsiBasin(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{
    spin_system_ptr = &spin_system;
    _max_counter = _get_max_counter(params);
    _init_E_objects();
    _init_S_objects();
}


_PsiBasinObjects* PsiBasin::_get_psi_basin_object_pointer(const std::string which)
{
    _PsiBasinObjects* psi_basin_object_pointer;
    if (which == "E") // This is the energy threshold ridge step
    {
        psi_basin_object_pointer = &E_objects;
    }
    else if (which == "S") // This is the entropic attractor ridge step
    {
        psi_basin_object_pointer = &S_objects;
    }
    else
    {
        throw std::runtime_error("Unknown 'which' (should be S, E)");
    }
    return psi_basin_object_pointer;
}


void PsiBasin::_step(const double current_waiting_time, const std::string which)
{

    _PsiBasinObjects* data_ptr = _get_psi_basin_object_pointer(which);
    if (!data_ptr->threshold_valid){return;}

    const parameters::StateProperties prev = spin_system_ptr->get_previous_state();
    const parameters::StateProperties curr = spin_system_ptr->get_current_state();

    const double prev_energy = prev.energy;
    const double curr_energy = curr.energy;

    // If the previous energy was above a threshold, we return immediately,
    // as we're either above a basin, or have just entered it
    if (prev_energy >= data_ptr->threshold){return;}

    // Append the current waiting time to the tracker
    data_ptr->_waiting_time += current_waiting_time;

    // Get the previous binary state representation as a string so that it's
    // compatible with our dictionaries
    const std::string state = spin_system_ptr->binary_state();

    // Insert this into the temporary configurations tracker- this is a set
    // which will not include duplicates
    data_ptr->_tmp_unique_configs_in_basin.insert(state);

    // Now we check the current energy. If it is >= than the threshold, we
    // need to perform some further logic, as this implies that the tracer
    // has just left a basin
    if (curr_energy < data_ptr->threshold){return;}

    // Get a copy of the current local waiting time
    const double local_waiting_time = data_ptr->_waiting_time;

    // Reset the waiting time
    data_ptr->_waiting_time = 0.0;

    // Get the key corresponding to the current local waiting time and append
    // that to the log scaled distribution of waiting times in basins
    const long long key = _get_key(local_waiting_time);
    if (key <= _max_counter - 1){data_ptr->_counter[key] += 1;}

    // Next, we get the number of unique configurations that have been in this
    // basin. This is simply the size of the set
    const long double n_unique = data_ptr->_tmp_unique_configs_in_basin.size();

    // This also needs to be converted into a key
    const long long key2 = _get_key(n_unique);

    // And again, this is now appended to its own counter
    if (key2 <= _max_counter - 1)
    {
        data_ptr->_counter_unique_configs_per_basin[key2] += 1;
    }

    // And finally, the temporary set is reset
    data_ptr->_tmp_unique_configs_in_basin.clear();
}

void PsiBasin::step(const double current_waiting_time)
{
    _step(current_waiting_time, "E");
    _step(current_waiting_time, "S");
}

void PsiBasin::_dump_outfile(const std::string which)
{
    _PsiBasinObjects* data_ptr = _get_psi_basin_object_pointer(which);
    if (!data_ptr->threshold_valid){return;}

    for (int ii=0; ii<_max_counter; ii++)
    {
        long long cc = data_ptr->_counter[ii];
        long long uc = data_ptr->_counter_unique_configs_per_basin[ii];
        fprintf(data_ptr->outfile, "%lli %lli\n", cc, uc);
    }
    fclose(data_ptr->outfile);
}

PsiBasin::~PsiBasin()
{
    _dump_outfile("E");
    _dump_outfile("S");
}
