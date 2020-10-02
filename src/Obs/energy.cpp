#include "Obs/energy.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

Energy::Energy(const std::string grid_directory)
{
    const std::string grid_location = grid_directory + "/energy.txt";
    load_long_long_grid_(grid, grid_location);
    length = grid.size();
    max_time = grid[length - 1];
}

void Energy::step_(const double simulation_clock, const Vals v)
{

    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}

    if (pointer > length - 1){return;}

    // Write to the outfile
    while (grid[pointer] < simulation_clock)
    {   
        fprintf(outfile, "%lli %lli %.05Lf %lli %.05Lf\n", grid[pointer],
            v.int_rep, v.energy, v.int_rep_IS, v.energy_IS);
        pointer += 1;
        if (pointer > length - 1){return;}
    }
}
