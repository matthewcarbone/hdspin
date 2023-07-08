#ifndef OBS1_H
#define OBS1_H

#include <vector>

#include "utils.h"
#include "spin.h"

class OnePointObservables
{
protected:
    const parameters::FileNames fnames;
    const parameters::SimulationParameters params;
    std::vector<long long> grid;
    int grid_length;
    const SpinSystem* spin_system_ptr;

    // Output files
    FILE* outfile_energy;
    FILE* outfile_energy_IS;
    FILE* outfile_capacity;
    FILE* outfile_acceptance_rate;
    FILE* outfile_inherent_structure_timings;
    FILE* outfile_walltime_per_waitingtime;

    // The pointer to the last-updated point on the grid
    unsigned int pointer = 0;

public:

    // Constructor: reads in the grid from the specified grid directory
    OnePointObservables(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);

    void step(const double simulation_clock);
    ~OnePointObservables();
};


#endif
