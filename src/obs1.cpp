#include "obs1.h"
#include "utils.h"

OnePointObservables::OnePointObservables(const parameters::FileNames fnames, parameters::SimulationParameters params, SpinSystem& spin_system) : fnames(fnames), params(params)
{
    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    grids::load_long_long_grid_(grid, grid_location);
    grid_length = grid.size();
    spin_system_ptr = &spin_system;

    // Energy
    outfile_energy = fopen(fnames.energy.c_str(), "w");

    // Cache capacity observable
    // First line is the total capacity
    const std::string cache_capacity_string = std::string(spin_system_ptr->get_emap_ptr()->get_capacity());
    outfile_capacity = fopen(fnames.cache_size.c_str(), "w");
    fprintf(outfile_capacity, "%s\n", cache_capacity_string.c_str());

    // Inherent structure calculation time observable
    outfile_acceptance_rate = fopen(fnames.acceptance_rate.c_str(), "w");

    // Acceptance rate observable
    outfile_inherent_structure_timings = fopen(fnames.inherent_structure_timings.c_str(), "w");
}

void OnePointObservables::step(const double simulation_clock)
{
    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}
    if (pointer > grid_length - 1){return;}

    // Energy
    const double energy = spin_system_ptr->get_previous_state().energy;

    // Get the current simulation statistics for some of the observables
    const parameters::SimulationStatistics sim_stats = spin_system_ptr->get_sim_stats();

    // Write to the outfile_capacity
    const std::string cache_size_string = std::string(spin_system_ptr->get_emap_ptr()->get_size());

    // Inherent structure calculation time
    const double denominiator = sim_stats.inherent_structure_calls ? sim_stats.inherent_structure_calls > 0 : 1;
    const double inherent_structure_time_per_timestep = sim_stats.inherent_structure_total_time / ((double) denominiator);

    // Get acceptance rates
    const double acceptance_rate = ((double) sim_stats.acceptances) / ((double) sim_stats.total_steps);

    while (grid[pointer] < simulation_clock)
    {   
        fprintf(outfile_energy, "%.08f\n", energy);
        fprintf(outfile_capacity, "%s\n", cache_size_string.c_str());
        fprintf(outfile_acceptance_rate, "%.08f\n", acceptance_rate);
        fprintf(outfile_inherent_structure_timings, "%.08f\n", inherent_structure_time_per_timestep);

        pointer += 1;
        if (pointer > grid_length - 1){return;}
    }
}

OnePointObservables::~OnePointObservables()
{
    fclose(outfile_energy);
    fclose(outfile_capacity);
    fclose(outfile_acceptance_rate);
    fclose(outfile_inherent_structure_timings);
}

