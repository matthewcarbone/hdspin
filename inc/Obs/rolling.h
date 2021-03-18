#ifndef OBS_ROLLING_H
#define OBS_ROLLING_H

#include "Obs/base.h"
#include "Utils/structures.h"

struct ridge_tracker
{
    long double mu0 = 0.0;
    long double mu1 = 0.0;
    long double S0 = 0.0;
    long double S1 = 0.0;
    long long counter = 1;
};

class Rolling : public Base
{
private:

    RuntimeParameters rtp;

    // Private data for the rolling means and variances for the ridge energy
    // calculation.
    ridge_tracker E_e_same, E_e_diff, S_e_same, S_e_diff;

    // Last energies that were under the threshold
    long double E_last_energy = 0.0;
    long double E_current_ridge = 0.0;
    long double S_last_energy = 0.0;
    long double S_current_ridge = 0.0;

    void _log_ridge_E(const Vals);
    void _log_ridge_S(const Vals);

public:

    // Constructor: reads in the grid from the specified grid directory
    Rolling(const FileNames, const RuntimeParameters);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step_(const Vals, const Vals);

    ~Rolling();
};


#endif
