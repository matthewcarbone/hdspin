#ifndef OBS_ENERGY_H
#define OBS_ENERGY_H

#include <vector>

#include "utils.h"
#include "spin.h"

class Energy
{
protected:

    parameters::FileNames fnames;
    parameters::SimulationParameters params;
    FILE *outfile;

    // The grid (in time) itself
    std::vector<long long> grid;

    // The length of the grid
    int grid_length;

    // The pointer to the last-updated point on the grid
    int pointer = 0;

    void _help_step(const long double simulation_clock, const double energy);

public:

    // Constructor: reads in the grid from the specified grid directory
    Energy(const parameters::FileNames fnames, const parameters::SimulationParameters params);

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param const long double simulation_clock The current timestep of the
     * simulation
     * @param const double energy The previous value for the energy (prev)
     */
    void step(const long double simulation_clock, const double energy);
    ~Energy();
};


#endif
