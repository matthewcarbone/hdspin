#ifndef OBS_CACHE_SIZE_H
#define OBS_CACHE_SIZE_H

#include <vector>

#include "utils.h"
#include "spin.h"

class CacheSize
{
protected:
    std::vector<long long> grid;
    int grid_length;
    FILE* outfile;
    EnergyMapping* emap_ptr;

    // The pointer to the last-updated point on the grid
    int pointer = 0;

public:

    // Constructor: reads in the grid from the specified grid directory
    CacheSize(const parameters::FileNames fnames, EnergyMapping& emap);

    void step(const long double simulation_clock);
    ~CacheSize();
};


#endif
