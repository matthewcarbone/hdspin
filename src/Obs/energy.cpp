#include "Obs/energy.h"
#include "Utils/utils.h"

EnergyBase::EnergyBase(const FileNames fnames,
    const RuntimeParameters rtp) : rtp(rtp)
{
    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    load_long_long_grid_(grid, grid_location);
    length = grid.size();
    max_time = grid[length - 1];
}

void EnergyBase::_help_step(const long double simulation_clock,
    const double energy)
{
    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}

    if (pointer > length - 1){return;}

    // Write to the outfile
    while (grid[pointer] < simulation_clock)
    {   
        fprintf(outfile, "%.05f\n", energy);

        pointer += 1;
        if (pointer > length - 1){return;}
    }
}



Energy::Energy(const FileNames fnames,
    const RuntimeParameters rtp) : EnergyBase(fnames, rtp)
{
    outfile = fopen(fnames.energy.c_str(), "w");
}

Energy::~Energy()
{
    fclose(outfile);
}

void Energy::step(const long double simulation_clock, const Vals v)
{
    _help_step(simulation_clock, v.energy);
}


EnergyAvgNeighbors::EnergyAvgNeighbors(const FileNames fnames,
    const RuntimeParameters rtp) : EnergyBase(fnames, rtp)
{
    outfile = fopen(fnames.energy_avg_neighbors.c_str(), "w");
}

EnergyAvgNeighbors::~EnergyAvgNeighbors()
{
    fclose(outfile);
}

void EnergyAvgNeighbors::step(
    const long double simulation_clock, const SpinSystem& sys)
{
    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}

    if (pointer > length - 1){return;}

    const double energy = sys.get_average_neighboring_energy();

    // Write to the outfile
    while (grid[pointer] < simulation_clock)
    {
        fprintf(outfile, "%.05f\n", energy);

        pointer += 1;
        if (pointer > length - 1){return;}
    }
}


EnergyInherentStructure::EnergyInherentStructure(const FileNames fnames,
    const RuntimeParameters rtp) : EnergyBase(fnames, rtp)
{
    if (rtp.memory != 0){outfile = fopen(fnames.energy_IS.c_str(), "w");}
}

EnergyInherentStructure::~EnergyInherentStructure()
{
    if (rtp.memory != 0){fclose(outfile);}
}

void EnergyInherentStructure::step(
    const long double simulation_clock, const Vals v)
{
    if (rtp.memory != 0)
    {
        _help_step(simulation_clock, v.energy_IS);
    }
}
