#ifndef OBS_ENERGY_H
#define OBS_ENERGY_H

#include <vector>

#include "Utils/structures.h"

class EnergyBase
{
protected:

    FileNames fnames;
    RuntimeParameters rtp;
    FILE *outfile;
    FILE *outfile_IS;

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
    EnergyBase(const FileNames, const RuntimeParameters);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void _help_step(const long double, const long long int_rep,
        const double energy);
};


class Energy : public EnergyBase
{
public:
    Energy(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals);
    ~Energy();
};


class EnergyInherentStructure : public EnergyBase
{
public:
    EnergyInherentStructure(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals);
    ~EnergyInherentStructure();
};


#endif
