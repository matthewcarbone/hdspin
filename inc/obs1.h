#ifndef OBS1_H
#define OBS1_H

#include <vector>

#include "utils.h"
#include "spin.h"

class OnePointObservables
{
protected:
    parameters::FileNames fnames;
    parameters::SimulationParameters params;
    std::vector<long long> grid;
    int grid_length;
    SpinSystem* spin_system_ptr;

    // Output files
    FILE* outfile_energy;
    FILE* outfile_capacity;
    FILE* outfile_acceptance_rate;
    FILE* outfile_inherent_structure_timings;

    // The pointer to the last-updated point on the grid
    int pointer = 0;

public:

    // Constructor: reads in the grid from the specified grid directory
    OnePointObservables(const parameters::FileNames fnames, parameters::SimulationParameters params, SpinSystem& spin_system);

    void step(const double simulation_clock);
    ~OnePointObservables();
};


#endif
