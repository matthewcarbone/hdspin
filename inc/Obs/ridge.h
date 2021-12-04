#ifndef OBS_RIDGE_H
#define OBS_RIDGE_H

#include "Utils/structures.h"

class RidgeBase
{
protected:

    FileNames fnames;
    RuntimeParameters rtp;

    FILE* outfile;

    // Last energies that were under the threshold
    double last_energy = 0.0;
    double current_ridge = 0.0;

    // Steps above the ridge energy
    int steps_above = 0;

    // Time above the ridge energy
    double time_above = 0.0;

    // The number of ridge energies logged
    int _logged = 0;

    // We must handle the case when the tracer STARTS above the threshold. In
    // this situation, we should not log the first time it drops below.
    bool _exited_first_basin = false;

    double _threshold;
    bool _threshold_valid = true;

    void _log_ridge(const double);

public:

    // Constructor: reads in the grid from the specified grid directory
    RidgeBase(const FileNames, const RuntimeParameters);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step_(const Vals, const Vals, const double);

    ~RidgeBase();
};



class RidgeEnergy : RidgeBase
{
public:
    RidgeEnergy(const FileNames, const RuntimeParameters);
};

class RidgeAttractor : RidgeBase
{
public:
    RidgeAttractor(const FileNames, const RuntimeParameters);
};


#endif
