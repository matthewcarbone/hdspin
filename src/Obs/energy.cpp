#include "Obs/energy.h"
#include "Utils/utils.h"

Energy::Energy(const FileNames fnames, const RuntimeParameters rtp) : rtp(rtp)
{
    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    load_long_long_grid_(grid, grid_location);
    length = grid.size();
    max_time = grid[length - 1];
    outfile = fopen(fnames.energy.c_str(), "w");
}

void Energy::step_(const long double simulation_clock, const Vals v)
{
    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}

    if (pointer > length - 1){return;}

    // Write to the outfile
    while (grid[pointer] < simulation_clock)
    {   
        if (rtp.memory != 0)
        {
            fprintf(outfile, "%lli %lli %.05f %lli %.05f\n", grid[pointer],
                v.int_rep, v.energy, v.int_rep_IS, v.energy_IS);
        }
        else
        {
            fprintf(outfile, "%lli %lli %.05f\n", grid[pointer], v.int_rep,
                v.energy);
        }

        pointer += 1;
        if (pointer > length - 1){return;}
    }
}

Energy::~Energy()
{
    fclose(outfile);
}
