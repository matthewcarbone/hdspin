/**
 * Two-point correlators such as psi and Pi.
 */

#include <math.h>

#include "obs2.h"
#include "utils.h"


size_t _get_max_counter(const utils::SimulationParameters params)
{
    // let N be the log base 10 total timesteps
    // log2(10^N) = N log2(10)
    // Give the max counter a lot of space
    return  params.log10_N_timesteps * log2(10.0) + 10;
}


/**
 * @brief Fills a long long vector with zeros up to _max_counter
 */
void _fill_counter(const size_t _max_counter, std::vector<long long> &counter)
{
    for (size_t ii=0; ii<_max_counter; ii++){counter.push_back(0);}
}


int _get_key(const double local_waiting_time)
{
    // If the waiting time is <= 1, round it to 1.
    if (local_waiting_time <= 1.0){return 0;} // 2^0

    // Otherwise return the log2-binned result
    // log2(1) == 0
    // log2(2) == 1 -> 1
    // log2(2.5) == 1.3 -> 2
    // log2(3) == 1.7 -> 2
    // log2(4) == 2 -> 2
    // log2(7) == 2.8 -> 3
    // log2(8) == 3 -> 3
    // log2(8.000001) == 3.000001 -> 4
    return ceil(log2(local_waiting_time));
}



// PSI CONFIG -----------------------------------------------------------------

PsiConfig::PsiConfig(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{
    spin_system_ptr = &spin_system;
    _max_counter = _get_max_counter(params);
    _fill_counter(_max_counter, _counter);
}

void PsiConfig::step(const double current_waiting_time)
{

    const utils::StateProperties prev = spin_system_ptr->get_previous_state();
    const utils::StateProperties curr = spin_system_ptr->get_current_state();

    // (x) ----------------------> ( )
    // prev prev prev prev prev   curr
    //
    //
    //
    //                              prev ...
    //                             (x)

    // No matter what, we update the internal states. If the configs are logged
    // to the counters they are reset in the helper functions.
    _waiting_time += current_waiting_time;

    // Now we attempt to log things. If the previous configuration representation
    // is different than the current one, we simply return.
    if (curr.state == prev.state){return;}

    // Otherwise we get the key
    const size_t key = _get_key(_waiting_time);

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


void PsiBasin::_init_E_data()
{
    _fill_counter(_max_counter, E_data._counter);
    _fill_counter(_max_counter, E_data._counter_unique_configs_per_basin);
    E_data.threshold = params.energetic_threshold;
    E_data.outfile = fopen(fnames.psi_basin_E.c_str(), "w");
}

void PsiBasin::_init_S_data()
{
    _fill_counter(_max_counter, S_data._counter);
    _fill_counter(_max_counter, S_data._counter_unique_configs_per_basin);
    if (!params.valid_entropic_attractor){S_data.threshold_valid = false;}
    S_data.threshold = params.entropic_attractor;
    if (S_data.threshold_valid)
    {
        S_data.outfile = fopen(fnames.psi_basin_S.c_str(), "w");
    }
}


PsiBasin::PsiBasin(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{
    spin_system_ptr = &spin_system;
    _max_counter = _get_max_counter(params);
    _init_E_data();
    _init_S_data();
}


PsiBasinData* PsiBasin::_get_psi_basin_data_pointer(const std::string which)
{
    if (which == "E") // This is the energy threshold
    {
        return &E_data;
    }
    else if (which == "S") // This is the entropic attractor
    {
        return &S_data;
    }
    else
    {
        throw std::runtime_error("Unknown 'which' (should be S, E)");
    }
}


void PsiBasin::_step(const double current_waiting_time, const std::string which)
{

    PsiBasinData* data_ptr = _get_psi_basin_data_pointer(which);
    if (!data_ptr->threshold_valid){return;}

    const utils::StateProperties prev = spin_system_ptr->get_previous_state();
    const utils::StateProperties curr = spin_system_ptr->get_current_state();

    const double prev_energy = prev.energy;
    const double curr_energy = curr.energy;

    // If the previous energy was above a threshold, we return immediately,
    // as we're either above a basin, or have just entered it
    if (prev_energy >= data_ptr->threshold){return;}

    // Append the current waiting time to the tracker
    data_ptr->_waiting_time += current_waiting_time;

    // Get the previous binary state representation as a string so that it's
    // compatible with our dictionaries
    const std::string state = spin_system_ptr->get_previous_state_string_rep();

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
    PsiBasinData* data_ptr = _get_psi_basin_data_pointer(which);
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



// PI CONFIG ------------------------------------------------------------------



Aging::Aging(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{
    spin_system_ptr = &spin_system;

    const std::string pi_1_grid_location = fnames.grids_directory + "/pi1.txt";
    utils::load_long_long_grid_(grid_pi1, pi_1_grid_location);
    length = grid_pi1.size();

    const std::string pi_2_grid_location = fnames.grids_directory + "/pi2.txt";
    utils::load_long_long_grid_(grid_pi2, pi_2_grid_location);
    const int grid_length2 = grid_pi2.size();

    assert(length == grid_length2);
}



// TODO - refactor!

void AgingConfig::_help_step_1(const double simulation_clock)
{

    // Two break conditions
    if (grid_pi1[pointer1] >= simulation_clock) {return;}
    if (pointer1 > length - 1) {return;}

    const std::string state = spin_system_ptr->get_previous_state_string_rep();

    while (grid_pi1[pointer1] < simulation_clock)
    {
        results1.push_back(state);
        pointer1 += 1;
        if (pointer1 > length - 1){break;}
    }
}


void AgingConfig::_help_step_2(const double simulation_clock)
{

    // Two break conditions
    if (grid_pi2[pointer2] >= simulation_clock) {return;}
    if (pointer2 > length - 1) {return;}

    const std::string state = spin_system_ptr->get_previous_state_string_rep();
    // std::cout << "state(2) " << state << std::endl;

    while (grid_pi2[pointer2] < simulation_clock)
    {
        results2.push_back(state);
        pointer2 += 1;
        if (pointer2 > length - 1){break;}
    }
}



AgingConfig::AgingConfig(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system) : Aging(fnames, params, spin_system)
{
    outfile = fopen(fnames.Pi_config.c_str(), "w");
}


void AgingConfig::step(const double simulation_clock)
{
    _help_step_1(simulation_clock);
    _help_step_2(simulation_clock);
}

AgingConfig::~AgingConfig()
{

    // Sanity checks to make sure that every vector was filled completely
    assert(length == results1.size());
    assert(length == results2.size());

    // Dump results to disk
    for (int ii = 0; ii < length; ii++)
    {
        const std::string s1 = results1[ii];
        const std::string s2 = results2[ii];
        fprintf(outfile, "%s %s\n", s1.c_str(), s2.c_str());
    }
    fclose(outfile);   
}


// PI BASIN -------------------------------------------------------------------


void AgingBasin::_init_E_data()
{
    E_data.threshold = params.energetic_threshold;
    E_data.outfile = fopen(fnames.Pi_basin_E.c_str(), "w");
}

void AgingBasin::_init_S_data()
{
    if (!params.valid_entropic_attractor){S_data.threshold_valid = false;}
    S_data.threshold = params.entropic_attractor;
    if (S_data.threshold_valid)
    {
        S_data.outfile = fopen(fnames.Pi_basin_S.c_str(), "w");
    }
}

AgingBasinData* AgingBasin::_get_psi_basin_data_pointer(const std::string which)
{
    if (which == "E") // This is the energy threshold
    {
        return &E_data;
    }
    else if (which == "S") // This is the entropic attractor
    {
        return &S_data;
    }
    else
    {
        throw std::runtime_error("Unknown 'which' (should be S, E)");
    }
}


AgingBasin::AgingBasin(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system) : Aging(fnames, params, spin_system)
{
    _init_E_data();
    _init_S_data();
}



void AgingBasin::_help_step_1_(const double simulation_clock, const double prev_energy, AgingBasinData* data_ptr)
{
    const int prev_state_in_basin = prev_energy < data_ptr->threshold ? 1 : 0;

    // Column 1: basin index
    // Column 2: in basin or not (1 if in basin, 0 if not in basin)
    while (grid_pi1[data_ptr->pointer1] < simulation_clock)
    {
        data_ptr->vec_basin_index_1.push_back(data_ptr->basin_index_1);
        data_ptr->vec_prev_state_in_basin_1.push_back(prev_state_in_basin);
        data_ptr->pointer1 += 1;
        if (data_ptr->pointer1 > length - 1){break;}
    }
}

void AgingBasin::_help_step_2_(const double simulation_clock, const double prev_energy, AgingBasinData* data_ptr)
{
    const int prev_state_in_basin = prev_energy < data_ptr->threshold ? 1 : 0;
    while (grid_pi2[data_ptr->pointer2] < simulation_clock)
    {
        data_ptr->vec_basin_index_2.push_back(data_ptr->basin_index_2);
        data_ptr->vec_prev_state_in_basin_2.push_back(prev_state_in_basin);
        data_ptr->pointer2 += 1;
        if (data_ptr->pointer2 > length - 1){break;}
    }
}


void AgingBasin::_help_step(const double simulation_clock, const std::string which)
{

    AgingBasinData* data_ptr = _get_psi_basin_data_pointer(which);
    if (!data_ptr->threshold_valid){return;}

    // Early return conditions
    if ((data_ptr->pointer1 > length - 1) && (data_ptr->pointer2 > length - 1)) {return;}

    // We'll now need the previous and current energies
    const utils::StateProperties prev = spin_system_ptr->get_previous_state();
    const utils::StateProperties curr = spin_system_ptr->get_current_state();

    const double prev_energy = prev.energy;
    const double curr_energy = curr.energy;

    if (data_ptr->pointer1 > length - 1){;}
    else
    {
        if (simulation_clock > grid_pi1[data_ptr->pointer1])
        {
            _help_step_1_(simulation_clock, prev_energy, data_ptr);
        }
        // Just entered a basin, iterate the basin index. We do this after
        // stepping the grids because the tracer is assumed to be in this
        // configuration until, but not including, the current time indexed by
        //the simulation clock.
        if ((prev_energy >= data_ptr->threshold) && (curr_energy < data_ptr->threshold))
        {
            data_ptr->basin_index_1++;
        }
    }

    if (data_ptr->pointer2 > length - 1){;}
    else
    {
        if (simulation_clock > grid_pi2[data_ptr->pointer2])
        {
            _help_step_2_(simulation_clock, prev_energy, data_ptr);
        }

        if ((prev_energy >= data_ptr->threshold) && (curr_energy < data_ptr->threshold))
        {
            data_ptr->basin_index_2++;
        }
    }
}

void AgingBasin::step(const double simulation_clock)
{
    _help_step(simulation_clock, "E");
    _help_step(simulation_clock, "S");
}


void AgingBasin::_dump_outfile(const std::string which)
{
    AgingBasinData* data_ptr = _get_psi_basin_data_pointer(which);
    if (!data_ptr->threshold_valid){return;}

    for (int ii=0; ii<length; ii++)
    {
        fprintf(data_ptr->outfile, "%lli %i %lli %i\n", data_ptr->vec_basin_index_1[ii],
            data_ptr->vec_prev_state_in_basin_1[ii], data_ptr->vec_basin_index_2[ii],
            data_ptr->vec_prev_state_in_basin_2[ii]);
    }
    fclose(data_ptr->outfile);
}

AgingBasin::~AgingBasin()
{
    _dump_outfile("E");
    _dump_outfile("S");
}
