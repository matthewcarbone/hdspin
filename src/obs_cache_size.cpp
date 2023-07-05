#include "obs_cache_size.h"
#include "utils.h"

CacheSize::CacheSize(const parameters::FileNames fnames, EnergyMapping& emap)
{
    const std::string grid_location = fnames.grids_directory + "/energy.txt";
    grids::load_long_long_grid_(grid, grid_location);
    grid_length = grid.size();
    outfile = fopen(fnames.cache_size.c_str(), "w");
    emap_ptr = &emap;

    // First line is the total capacity
    const std::string cache_capacity_string = std::string(emap_ptr->get_capacity());
    fprintf(outfile, "%s\n", cache_capacity_string.c_str());
}

void CacheSize::step(const long double simulation_clock)
{
    // No updates necessary
    if (simulation_clock <= grid[pointer]){return;}

    if (pointer > grid_length - 1){return;}

    // Write to the outfile
    const std::string cache_size_string = std::string(emap_ptr->get_size());
    while (grid[pointer] < simulation_clock)
    {   
        fprintf(outfile, "%s\n", cache_size_string.c_str());

        pointer += 1;
        if (pointer > grid_length - 1){return;}
    }
}

CacheSize::~CacheSize()
{
    fclose(outfile);
}

