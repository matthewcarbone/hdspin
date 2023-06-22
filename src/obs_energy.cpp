#include "obs_energy.h"
#include "utils.h"

Energy::Energy(const parameters::FileNames fnames,
    const parameters::SimulationParameters params) : params(params)
{
    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    grids::load_long_long_grid_(grid, grid_location);
    grid_length = grid.size();
    outfile = fopen(fnames.energy.c_str(), "w");
}

void Energy::_help_step(const long double simulation_clock, const double energy)
{
    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}

    if (pointer > grid_length - 1){return;}

    // Write to the outfile
    while (grid[pointer] < simulation_clock)
    {   
        fprintf(outfile, "%.05f\n", energy);

        pointer += 1;
        if (pointer > grid_length - 1){return;}
    }
}

Energy::~Energy()
{
    fclose(outfile);
}

void Energy::step(const long double simulation_clock, const double energy)
{
    _help_step(simulation_clock, energy);
}
