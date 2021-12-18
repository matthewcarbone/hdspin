#ifndef OBS_RIDGE_H
#define OBS_RIDGE_H

#include <unordered_set>

#include "Utils/structures.h"

class RidgeBase
{
protected:

    FileNames fnames;
    RuntimeParameters rtp;

    FILE* outfile;

    // Last energies that were under the threshold
    double _last_energy = 0.0;
    double _current_ridge = 0.0;

    // Steps above the ridge energy
    int _steps_above = 0;
    std::unordered_set<int> _unique_configs_above;

    // Time above the ridge energy
    double _time_above = 0.0;

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
    void step(const Vals, const Vals, const double);

    ~RidgeBase();
};



class RidgeEnergy : public RidgeBase
{
public:
    RidgeEnergy(const FileNames, const RuntimeParameters);
};

class RidgeAttractor : public RidgeBase
{
public:
    RidgeAttractor(const FileNames, const RuntimeParameters);
};


#endif
