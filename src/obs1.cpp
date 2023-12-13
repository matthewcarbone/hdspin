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

OnePointObservables::OnePointObservables(const utils::SimulationParameters params, const SpinSystem& spin_system) : params(params)
{

    spin_system_ptr = &spin_system;

    utils::load_long_long_grid_(grid, ENERGY_GRID_PATH);
    grid_length = grid.size();

    // Cache capacity observable
    // First line is the total capacity
    const std::string cache_capacity_string = std::string(spin_system_ptr->get_emap_ptr()->get_capacity());

    // DON't FORGET TO WRITE THE CACHE CAPACITY
    //fprintf(outfile_capacity, "%s\n", cache_capacity_string.c_str());

    // Initialize the ridge energies from the provided parameters
    ridge_E_object.threshold = params.energetic_threshold;
    ridge_E_object.threshold_valid = true;
    ridge_S_object.threshold = params.entropic_attractor;
    ridge_S_object.threshold_valid = params.valid_entropic_attractor;
}


RidgeEnergyObject* OnePointObservables::_get_RidgeEnergyObject_ptr(const std::string which_ridge)
{
    RidgeEnergyObject* ridge_ptr;
    if (which_ridge == "E") // This is the energy threshold ridge step
    {
        ridge_ptr = &ridge_E_object;
    }
    else if (which_ridge == "S") // This is the entropic attractor ridge step
    {
        ridge_ptr = &ridge_S_object;
    }
    else
    {
        throw std::runtime_error("Unknown ridge threshold");
    }
    return ridge_ptr;
}


void OnePointObservables::_step_ridge(const double waiting_time, const double simulation_clock, const std::string which_ridge)
{

    RidgeEnergyObject* ridge_ptr = _get_RidgeEnergyObject_ptr(which_ridge);

    if (!ridge_ptr->threshold_valid){return;}

    const utils::StateProperties prev = spin_system_ptr->get_previous_state();
    const utils::StateProperties curr = spin_system_ptr->get_current_state();

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

void OnePointObservables::_log_ridge(const std::string which_ridge)
{

    RidgeEnergyObject* ridge_ptr = _get_RidgeEnergyObject_ptr(which_ridge);

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

    const double median = ridge_ptr->streaming_median.median();
    ridge_ptr->vec_medians.push_back(median);

    const double mean = ridge_ptr->streaming_mean.mean();
    ridge_ptr->vec_means.push_back(mean);

    ridge_ptr->vec_total_steps.push_back(ridge_ptr->total_steps);
}

void OnePointObservables::step(const double waiting_time, const double simulation_clock)
{

    // No matter what we step the ridge energies
    _step_ridge(waiting_time, simulation_clock, "E");
    _step_ridge(waiting_time, simulation_clock, "S");

    // And then continue on...
    const utils::StateProperties prev = spin_system_ptr->get_previous_state(); 

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
    const utils::SimulationStatistics sim_stats = spin_system_ptr->get_sim_stats();

    // Write to the outfile_capacity
    const std::string cache_size_string = std::string(spin_system_ptr->get_emap_ptr()->get_size());

    // Get acceptance rates
    const double acceptance_rate = ((double) sim_stats.acceptances) / ((double) sim_stats.total_steps);

    while (grid[pointer] < simulation_clock)
    {   
        vec_energy.push_back(energy);
        
        if (params.calculate_inherent_structure_observables)
        {
            vec_energy_IS.push_back(energy_IS);
        }
        vec_cache_size.push_back(cache_size_string.c_str());
        vec_acceptance_rate.push_back(acceptance_rate);
        vec_walltime_per_waiting_time.push_back(
            sim_stats.total_wall_time/sim_stats.total_waiting_time);
        _log_ridge("E");
        _log_ridge("S");

        pointer += 1;
        if (pointer > grid_length - 1){return;}
    }
}

json OnePointObservables::as_json() const {
    json j;
    j["energy"] = vec_energy;
    j["energy_inherent_structure"] = vec_energy_IS;
    j["cache_size"] = vec_cache_size;
    j["acceptance_rate"] = vec_acceptance_rate;
    j["walltime_per_waitingtime"] = vec_walltime_per_waiting_time;
    j["ridge_E_mean"] = ridge_E_object.vec_means;
    j["ridge_E_median"] = ridge_E_object.vec_medians;
    j["ridge_E_total_steps"] = ridge_E_object.vec_total_steps;
    j["ridge_S_mean"] = ridge_S_object.vec_means;
    j["ridge_S_median"] = ridge_S_object.vec_medians;
    j["ridge_S_total_steps"] = ridge_S_object.vec_total_steps;
    return j;
}
