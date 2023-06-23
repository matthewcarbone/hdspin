#ifndef OBS_PSI_H
#define OBS_PSI_H

#include <vector>

#include "Utils/structures.h"

class PsiConfigBase
{
protected:

    FileNames fnames;
    FILE *outfile;

    // Runtime parameters
    RuntimeParameters rtp;

    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    long long _max_counter;

    std::vector<long long> _counter;

    // Keep track internally of the waiting time for both the standard
    // trajectory and the inherent structure
    long double _waiting_time = 0.0;

    // Helper methods
    void _help_step();

public:

    PsiConfigBase(const FileNames, const RuntimeParameters);
};



class PsiConfig : public PsiConfigBase
{
public:
    PsiConfig(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
    ~PsiConfig();
};


class PsiConfigInherentStructure : public PsiConfigBase
{
public:
    PsiConfigInherentStructure(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
    ~PsiConfigInherentStructure();
};




class PsiBasinBase
{
protected:

    FileNames fnames;
    FILE *outfile;

    bool inherent_structure, energetic_barrier;

    // Runtime parameters
    RuntimeParameters rtp;

    // The threshold energy in question
    double _threshold;
    bool _threshold_valid = true;

    std::vector<long long> _counter;

    // Number of unique configs per basin
    std::vector<long long> _counter_unique_configs_per_basin;
    std::vector<long long> _tmp_unique_configs_in_basin;

    long double _waiting_time = 0.0;

    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    long long _max_counter;

    void _help_step(const double, const double, const long long,
        const long double);
    void _dump_outfile();

public:
    PsiBasinBase(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
};


class PsiBasinThreshold : public PsiBasinBase
{
public:
    PsiBasinThreshold(const FileNames, const RuntimeParameters);
    ~PsiBasinThreshold();
};


class PsiBasinAttractor : public PsiBasinBase
{
public:
    PsiBasinAttractor(const FileNames, const RuntimeParameters);
    ~PsiBasinAttractor();
};


class PsiBasinThresholdInherentStructure : public PsiBasinBase
{
public:
    PsiBasinThresholdInherentStructure(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
    ~PsiBasinThresholdInherentStructure();
};


class PsiBasinAttractorInherentStructure : public PsiBasinBase
{
public:
    PsiBasinAttractorInherentStructure(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
    ~PsiBasinAttractorInherentStructure();
};


#endif
