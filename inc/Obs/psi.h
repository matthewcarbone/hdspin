#ifndef OBS_PSI_H
#define OBS_PSI_H

#include "Obs/base.h"
#include "Utils/structures.h"

class PsiConfig : public Base
{
private:

    // Runtime parameters
    RuntimeParameters rtp;

    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    long long max_counter;

    // The grid (in time) itself
    std::vector<long long> counter;
    std::vector<long long> counter_IS;

    // Keep track internally of the waiting time for both the standard
    // trajectory and the inherent structure
    long double waiting_time = 0.0;
    long double waiting_time_IS = 0.0;

    // Helper methods
    void _help_step(const bool);

public:

    PsiConfig(const FileNames, const RuntimeParameters);
    void step_(const long double, const Vals, const Vals);
    ~PsiConfig();
};


#endif
