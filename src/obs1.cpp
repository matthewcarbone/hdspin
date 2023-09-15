#include "obs1.h"
#include "utils.h"


StreamingMedian::StreamingMedian(){}

void StreamingMedian::update(const double v)
{
    if(max_heap.empty() || max_heap.top() >= v) max_heap.push(v);
    else min_heap.push(v);

    if(max_heap.size() > min_heap.size() + 1){
        min_heap.push(max_heap.top());
        max_heap.pop();
    }
    else if(max_heap.size() < min_heap.size()){
        max_heap.push(min_heap.top());
        min_heap.pop();
    }
}

double StreamingMedian::median() const
{
    if (max_heap.size() == 0 & min_heap.size() == 0)
    {
        return 0.0;
    }
    if(max_heap.size() == min_heap.size()){
        return max_heap.top()/2.0 + min_heap.top()/2.0;
    }
    return max_heap.top();
}

OnePointObservables::OnePointObservables(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{

    spin_system_ptr = &spin_system;

    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    grids::load_long_long_grid_(grid, grid_location);
    grid_length = grid.size();
    spin_system_ptr = &spin_system;

    // Energy
    outfile_energy = fopen(fnames.energy.c_str(), "w");

    // Inherent structure energy
    if (params.calculate_inherent_structure_observables)
    {
        outfile_energy_IS = fopen(fnames.energy_IS.c_str(), "w");
    }
    
    // Cache capacity observable
    // First line is the total capacity
    const std::string cache_capacity_string = std::string(spin_system_ptr->get_emap_ptr()->get_capacity());
    outfile_capacity = fopen(fnames.cache_size.c_str(), "w");
    fprintf(outfile_capacity, "%s\n", cache_capacity_string.c_str());

    // Acceptance rate
    outfile_acceptance_rate = fopen(fnames.acceptance_rate.c_str(), "w");

    // Wall time/timestep
    outfile_walltime_per_waitingtime = fopen(fnames.walltime_per_waitingtime.c_str(), "w");

    // Ridge energy
    ridge_E_objects.outfile = fopen(fnames.ridge_E.c_str(), "w");
    ridge_E_objects.threshold = params.energetic_threshold;
    ridge_S_objects.threshold_valid = params.valid_entropic_attractor;
    ridge_S_objects.threshold = params.entropic_attractor;
    if (ridge_S_objects.threshold_valid)
    {
        ridge_S_objects.outfile = fopen(fnames.ridge_S.c_str(), "w");
    }
}


_RidgeEnergyObjects* OnePointObservables::_get_ridge_pointer(const std::string which_ridge)
{
    _RidgeEnergyObjects* ridge_ptr;
    if (which_ridge == "E") // This is the energy threshold ridge step
    {
        ridge_ptr = &ridge_E_objects;
    }
    else if (which_ridge == "S") // This is the entropic attractor ridge step
    {
        ridge_ptr = &ridge_S_objects;
    }
    else
    {
        throw std::runtime_error("Unknown ridge threshold");
    }
    return ridge_ptr;
}


void OnePointObservables::_step_ridge(const double waiting_time, const double simulation_clock, const std::string which_ridge)
{

    _RidgeEnergyObjects* ridge_ptr = _get_ridge_pointer(which_ridge);

    if (!ridge_ptr->threshold_valid){return;}

    const parameters::StateProperties prev = spin_system_ptr->get_previous_state();
    const parameters::StateProperties curr = spin_system_ptr->get_current_state();

    const double _prev_energy = prev.energy;
    const double _curr_energy = curr.energy;

    // Just exited a basin, track the previous energy
    if ((prev.energy < ridge_ptr->threshold) && (curr.energy >= ridge_ptr->threshold))
    {
        ridge_ptr->last_energy = prev.energy;
        ridge_ptr->current_ridge = _curr_energy;
        ridge_ptr->exited_first_basin = true;
    }

    // Still above the basin. The current ridge energy is defined as the
    // maximum energy reached above the threshold.
    else if ((prev.energy >= ridge_ptr->threshold) && (curr.energy >= ridge_ptr->threshold))
    {
        ridge_ptr->current_ridge = _curr_energy > ridge_ptr->current_ridge ? _curr_energy : ridge_ptr->current_ridge;
    }

    // Just dropped back below the threshold: log the current ridge energy
    else if ((prev.energy >= ridge_ptr->threshold) && (curr.energy < ridge_ptr->threshold))
    {
        if (ridge_ptr->exited_first_basin)
        {
            ridge_ptr->streaming_median.update(ridge_ptr->current_ridge);
            ridge_ptr->streaming_mean.update(ridge_ptr->current_ridge);
            ridge_ptr->total_steps += 1;
        }
    }
}

void OnePointObservables::_ridge_writeout(const std::string which_ridge)
{

    _RidgeEnergyObjects* ridge_ptr = _get_ridge_pointer(which_ridge);

    if (!ridge_ptr->threshold_valid){return;}

    double total_steps;
    if (ridge_ptr->total_steps > 0)
    {
        total_steps = ((double) ridge_ptr->total_steps);
    }
    else
    {
        total_steps = 1;
    }
    const double _median = ridge_ptr->streaming_median.median();
    const double _mean = ridge_ptr->streaming_mean.mean();
    fprintf(ridge_ptr->outfile, "%.08f %.08f %lli\n", _mean, _median, ridge_ptr->total_steps);
}

void OnePointObservables::step(const double waiting_time, const double simulation_clock)
{

    // No matter what we step the ridge energies
    _step_ridge(waiting_time, simulation_clock, "E");
    _step_ridge(waiting_time, simulation_clock, "S");

    // And then continue on...
    const parameters::StateProperties prev = spin_system_ptr->get_previous_state(); 

    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}
    if (pointer > grid_length - 1){return;}

    // Energy
    const double energy = prev.energy;

    // Energy inherent structure
    double energy_IS;
    if (params.calculate_inherent_structure_observables)
    {
        const ap_uint<PRECISON> inherent_structure = spin_system_ptr->get_emap_ptr()->get_inherent_structure(prev.state);
        energy_IS = spin_system_ptr->get_emap_ptr()->get_config_energy(inherent_structure);
    }
    else
    {
        energy_IS = NAN;
    }
    

    // Get the current simulation statistics for some of the observables
    const parameters::SimulationStatistics sim_stats = spin_system_ptr->get_sim_stats();

    // Write to the outfile_capacity
    const std::string cache_size_string = std::string(spin_system_ptr->get_emap_ptr()->get_size());

    // Get acceptance rates
    const double acceptance_rate = ((double) sim_stats.acceptances) / ((double) sim_stats.total_steps);

    while (grid[pointer] < simulation_clock)
    {   
        fprintf(outfile_energy, "%.08f\n", energy);
        if (params.calculate_inherent_structure_observables)
        {
            fprintf(outfile_energy_IS, "%.08f\n", energy_IS);
        }
        fprintf(outfile_capacity, "%s\n", cache_size_string.c_str());
        fprintf(outfile_acceptance_rate, "%.08f\n", acceptance_rate);
        fprintf(outfile_walltime_per_waitingtime, "%.08f\n", sim_stats.total_wall_time/sim_stats.total_waiting_time);
        _ridge_writeout("E");
        _ridge_writeout("S");

        pointer += 1;
        if (pointer > grid_length - 1){return;}
    }
}

OnePointObservables::~OnePointObservables()
{
    fclose(outfile_energy);
    if (params.calculate_inherent_structure_observables)
    {
        fclose(outfile_energy_IS);
    }
    fclose(outfile_capacity);
    fclose(outfile_acceptance_rate);
    fclose(outfile_walltime_per_waitingtime);
    if (ridge_S_objects.threshold_valid)
    {
        fclose(ridge_S_objects.outfile);
    }
}

