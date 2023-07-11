#include "obs1.h"
#include "utils.h"


ObsBase::ObsBase(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{
    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    grids::load_long_long_grid_(grid, grid_location);
    grid_length = grid.size();
    spin_system_ptr = &spin_system;
};

RidgeBase::RidgeBase(const parameters::FileNames fnames,
    const parameters::SimulationParameters params, const SpinSystem& spin_system) : ObsBase(fnames, params, spin_system){}

void RidgeBase::step(const double waiting_time, const double simulation_clock)
{

    if (!_threshold_valid){return;}

    const parameters::StateProperties prev = spin_system_ptr->get_previous_state();
    const parameters::StateProperties curr = spin_system_ptr->get_current_state();

    const double _prev_energy = prev.energy;
    const double _curr_energy = curr.energy;

    // std::cout << prev.energy << " " << _threshold << std::endl;

    // Just exited a basin, track the previous energy
    if ((prev.energy < _threshold) && (curr.energy >= _threshold))
    {
        // The energy before exiting the basin is defined as the inherent
        // structure energy if we are using a system with memory. Else, we just
        // use the energy itself
        // _last_energy = (rtp.memory != 0) ? prev.energy_IS : prev.energy
        _last_energy = prev.energy;
        _current_ridge = _curr_energy;
        // const std::string str_rep = state::string_rep_from_arbitrary_precision_integer(curr.state, params.N_spins);
        // _unique_configs_above.insert(str_rep);
        _exited_first_basin = true;
    }

    // Still above the basin. The current ridge energy is defined as the
    // maximum energy reached above the threshold.
    else if ((prev.energy >= _threshold) && (curr.energy >= _threshold))
    {

        _current_ridge = _curr_energy > _current_ridge ? _curr_energy : _current_ridge;
        // _current_ridge = curr.energy > _current_ridge ? curr.energy : _current_ridge;
        // _steps_above += 1;
        // const std::string str_rep = state::string_rep_from_arbitrary_precision_integer(curr.state, params.N_spins);
        // _unique_configs_above.insert(str_rep);
        // _time_above += waiting_time;
    }

    // Just dropped back below the threshold: log the current ridge energy
    else if ((prev.energy >= _threshold) && (curr.energy < _threshold))
    {
        // std::cout << _exited_first_basin << std::endl;
        // _steps_above += 1;
        // _time_above += waiting_time;
        if (_exited_first_basin)
        {
            // const double curr.energy_compare = curr.energy;
            // _log_ridge(curr.energy_compare, simulation_clock);
            // std::cout << prev.energy << " " <<  curr.energy << std::endl;
            _ridge_energy_accumulator += _current_ridge;
            _total_steps += 1;
        }
        // _steps_above = 0;
        // _time_above = 0.0;
        // _unique_configs_above.clear();
    }

    if (grid[pointer] < simulation_clock)
    {
        double total_steps;
        if (_total_steps > 0)
        {
            total_steps = ((double) _total_steps);
        }
        else
        {
            total_steps = 1;
        }
        const double average = _ridge_energy_accumulator / total_steps;
        while (grid[pointer] < simulation_clock)
        {   
            fprintf(outfile, "%.08f %lli\n", average, _total_steps);
            pointer += 1;
            if (pointer > grid_length - 1){return;}
        }
    }
}

RidgeBase::~RidgeBase()
{
    if (_threshold_valid){fclose(outfile);}
}


RidgeE::RidgeE(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : RidgeBase(fnames, params, spin_system)
{
    outfile = fopen(fnames.ridge_E.c_str(), "w");
    _threshold = params.energetic_threshold;
}


RidgeS::RidgeS(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : RidgeBase(fnames, params, spin_system)
{
    _threshold_valid = params.valid_entropic_attractor;
    _threshold = params.entropic_attractor;
    if (_threshold_valid)
    {
        outfile = fopen(fnames.ridge_S.c_str(), "w");
    }
}

OnePointObservables::OnePointObservables(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : ObsBase(fnames, params, spin_system)
{
    // Energy
    outfile_energy = fopen(fnames.energy.c_str(), "w");

    // Inherent structure energy
    outfile_energy_IS = fopen(fnames.energy_IS.c_str(), "w");

    // Cache capacity observable
    // First line is the total capacity
    const std::string cache_capacity_string = std::string(spin_system_ptr->get_emap_ptr()->get_capacity());
    outfile_capacity = fopen(fnames.cache_size.c_str(), "w");
    fprintf(outfile_capacity, "%s\n", cache_capacity_string.c_str());

    // Inherent structure calculation time observable
    outfile_acceptance_rate = fopen(fnames.acceptance_rate.c_str(), "w");

    // Wall time/timestep
    outfile_walltime_per_waitingtime = fopen(fnames.walltime_per_waitingtime.c_str(), "w");
}

void OnePointObservables::step(const double waiting_time, const double simulation_clock)
{

    const parameters::StateProperties prev = spin_system_ptr->get_previous_state(); 

    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}
    if (pointer > grid_length - 1){return;}

    // Energy
    const double energy = prev.energy;

    // Energy inherent structure
    const ap_uint<PRECISON> inherent_structure = spin_system_ptr->get_emap_ptr()->get_inherent_structure(prev.state);
    const double energy_IS = spin_system_ptr->get_emap_ptr()->get_config_energy(inherent_structure);

    // Get the current simulation statistics for some of the observables
    const parameters::SimulationStatistics sim_stats = spin_system_ptr->get_sim_stats();

    // Write to the outfile_capacity
    const std::string cache_size_string = std::string(spin_system_ptr->get_emap_ptr()->get_size());

    // Get acceptance rates
    const double acceptance_rate = ((double) sim_stats.acceptances) / ((double) sim_stats.total_steps);

    while (grid[pointer] < simulation_clock)
    {   
        fprintf(outfile_energy, "%.08f\n", energy);
        fprintf(outfile_energy_IS, "%.08f\n", energy_IS);
        fprintf(outfile_capacity, "%s\n", cache_size_string.c_str());
        fprintf(outfile_acceptance_rate, "%.08f\n", acceptance_rate);
        fprintf(outfile_walltime_per_waitingtime, "%.08f\n", sim_stats.total_wall_time/sim_stats.total_waiting_time);

        pointer += 1;
        if (pointer > grid_length - 1){return;}
    }
}

OnePointObservables::~OnePointObservables()
{
    fclose(outfile_energy);
    fclose(outfile_energy_IS);
    fclose(outfile_capacity);
    fclose(outfile_acceptance_rate);
    fclose(outfile_walltime_per_waitingtime);
}

