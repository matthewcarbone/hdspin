#include <iostream>
#include <vector>
#include <fstream>

#include "general_utils.h"
#include "grid_utils.h"


// ============================================================================
// Energy =====================================================================
// ============================================================================

EnergyGrid::EnergyGrid(const std::string grid_location)
{
    // printf("Initializing grid from %s\n", grid_location.c_str());
    std::ifstream myfile (grid_location);
    std::string line;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            grid.push_back(stoll(line));
        }
        myfile.close();
    }
    length = grid.size();
    max_time = grid[length - 1];
    // printf("Grid initialized with min/max times %lli/%lli and len %i\n",
    //     grid[0], max_time, length);
}

std::vector<long long> EnergyGrid::get_grid(){return grid;}

void EnergyGrid::open_outfile(const std::string d)
{
    outfile = fopen(d.c_str(), "w");
}

void EnergyGrid::close_outfile(){fclose(outfile);}

void EnergyGrid::step(const double current_time, const double current_energy,
    const int *config, const int N_spins, const double *energy_array,
    long long *inherent_structure_mapping)
{

    // No updates necessary
    if (current_time <= grid[pointer]){return;}

    if (pointer > length - 1){return;}

    // Get the current configuration integer representations and energies
    const long long config_int = binary_vector_to_int(config, N_spins);
    
    // Update the inherent structure dictionary
    long long config_IS_int;
    if (inherent_structure_mapping[config_int] != -1)
    {
        // We've computed the inherent structure before, no need to do
        // it again
        config_IS_int = inherent_structure_mapping[config_int];
    }
    else
    {
        config_IS_int = compute_inherent_structure(config, energy_array,
            N_spins);
        inherent_structure_mapping[config_int] = config_IS_int;
    }
    const double energy_IS = energy_array[config_IS_int];

    // Write to the outfile
    while (grid[pointer] < current_time)
    {   
        fprintf(outfile, "%lli %lli %.05f %lli %.05f\n", grid[pointer],
            config_int, current_energy, config_IS_int, energy_IS);
        pointer += 1;
        if (pointer > length - 1){return;}
    }
}
