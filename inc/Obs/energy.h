#ifndef OBS_ENERGY_H
#define OBS_ENERGY_H

#include "Obs/base.h"
#include "Utils/structures.h"

class Energy : public Base
{
private:

    // The maximum time on the grid
    long long max_time;

    // The grid (in time) itself
    std::vector<long long> grid;

    // The length of the grid
    int length;

    // The pointer to the last-updated point on the grid
    int pointer = 0;

public:

    // Constructor: reads in the grid from the specified grid directory
    Energy(const FileNames);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step_(const long double, const Vals);
};


#endif
