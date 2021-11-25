#ifndef OBS_RIDGE_H
#define OBS_RIDGE_H

#include "Obs/base.h"
#include "Utils/structures.h"

struct ridge_tracker
{
    double mu0 = 0.0;
    double mu1 = 0.0;
    double S0 = 0.0;
    double S1 = 0.0;
    long long counter = 1;
    double current_min = 1e15;
    double current_max = -1e15;
};

class RidgeEnergy : public Base
{
private:

    // inherent_structure
    // 0 if using all standard trajectories
    // 1 if using all IS trajectories
    // 2 if using standard for ridge calculation but IS for diff/same
    int inherent_structure;

    // energetic_threshold
    // True or false for E/S threshold usage.
    bool energetic_threshold;

    bool log_all;

    double threshold;

    RuntimeParameters rtp;

    // Private data for the rolling means and variances for the ridge energy
    // calculation.
    ridge_tracker same, diff;

    // Last energies that were under the threshold
    double last_energy = 0.0;
    double current_ridge = 0.0;

    // Steps above the ridge energy
    int steps_above = 0;

    // Time above the ridge energy
    double time_above = 0.0;

    // We must handle the case when the tracer STARTS above the threshold. In
    // this situation, we should not log the first time it drops below.
    bool exited_first_basin = false;

    void _log_ridge_(const double);

public:

    // Constructor: reads in the grid from the specified grid directory
    RidgeEnergy(const FileNames, const RuntimeParameters, const int,
        const bool, const bool);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step_(const Vals, const Vals, const double);

    ~RidgeEnergy();
};


#endif
