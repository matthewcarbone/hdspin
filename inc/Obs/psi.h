#ifndef OBS_PSI_H
#define OBS_PSI_H

#include "Obs/base.h"
#include "Utils/structures.h"

class PsiConfig : public Base
{
private:

    bool inherent_structure;

    // Runtime parameters
    RuntimeParameters rtp;

    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    long long max_counter;

    std::vector<long long> counter;

    // Keep track internally of the waiting time for both the standard
    // trajectory and the inherent structure
    long double waiting_time = 0.0;

    // Helper methods
    void _help_step_();

public:

    PsiConfig(const FileNames, const RuntimeParameters, const bool);
    void step_(const long double, const Vals, const Vals);
    ~PsiConfig();
};


class PsiBasin : public Base
{
private:

    bool inherent_structure, energetic_barrier;

    // Runtime parameters
    RuntimeParameters rtp;

    // The threshold energy in question
    double threshold;

    std::vector<long long> counter;

    // Number of unique configs per basin
    std::vector<long long> counter_unique_configs_per_basin;
    std::vector<long long> tmp_unique_configs_in_basin;

    long double waiting_time = 0.0;

    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    long long max_counter;

public:
    PsiBasin(const FileNames, const RuntimeParameters, const bool, const bool);
    void step_(const long double, const Vals, const Vals);
    ~PsiBasin();
};


#endif
