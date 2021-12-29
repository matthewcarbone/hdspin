#include "Obs/configs.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

Configs::Configs(const FileNames fnames,
    const RuntimeParameters rtp) : rtp(rtp), fnames(fnames)
{
    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    load_long_long_grid_(grid, grid_location);
    length = grid.size();
    max_time = grid[length - 1];
    outfile = fopen(fnames.unique_configs.c_str(), "w");
}

void Configs::step(const long double simulation_clock,
    const Vals prev)
{

    // Update the unique configs unordered set
    _unique_configs.insert(prev.int_rep);

    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}
    if (pointer > length - 1){return;}

    // Write to the outfile
    while (grid[pointer] < simulation_clock)
    {
        int size = _unique_configs.size();
        fprintf(outfile, "%i\n", size);

        pointer += 1;
        if (pointer > length - 1){return;}
    }
}

Configs::~Configs()
{
    fclose(outfile);
}
