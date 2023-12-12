/**
 * Max energy over ranges t/2 -> t
 */

#include <math.h>

#include "obs2.h"
#include "utils.h"
#include "emax.h"


EMaxt2::EMaxt2(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system) : fnames(fnames), params(params)
{
    spin_system_ptr = &spin_system;

    // We need the pi2 grid since this goes from 0 to the length of the
    // simulation
    std::vector<long long> grid2;
    const std::string grid_location = fnames.grids_directory + "/pi2.txt";
    utils::load_long_long_grid_(grid2, grid_location);
    length = grid2.size();

    // Grid 1 is grid2 // 2
    for (size_t ii=0; ii<length; ii++)
    {
        EMaxt2Data d;
        d.tmax = grid2[ii];
        d.tmin = d.tmax / 2;
        trackers[ii] = d;
    }
}

void EMaxt2::step(const double simulation_clock) const
{
    const utils::StateProperties prev = spin_system_ptr->get_previous_state();
    const double energy = prev.energy;

    std::vector<size_t> pop_these;
    for (std::pair<const size_t, EMaxt2Data>& n : trackers)
    {
        size_t key = n.first;
        EMaxt2Data* d = &n.second;

        // If the simulation clock is greater than the minimum time, we compare
        // the current energy value to the current max energy
        if (simulation_clock >= d->tmin)
        {
            d->max_energy = fmax(energy, d->max_energy);    
        }
        
        // If the simulation clock is greater than the maximum time, we can
        // actually pop it from the trackers after recording its max energy
        if (simulation_clock >= d->tmax)
        {
            pop_these.push_back(key);
        }
    }

    for (size_t ii=0; ii<pop_these.size(); ii++)
    {
        const int key = pop_these[ii];
        max_energies[key] = trackers[key].max_energy;
        trackers.erase(key);
    }
}


EMaxt2::~EMaxt2()
{
    outfile = fopen(fnames.max_energy.c_str(), "w");
    for (size_t ii=0; ii<length; ii++)
    {
        fprintf(outfile, "%.08f\n", max_energies[ii]);
    }
    fclose(outfile);
}
