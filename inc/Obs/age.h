#ifndef OBS_AGE_H
#define OBS_AGE_H

#include <vector>

#include "Utils/structures.h"


class Aging
{
protected:

    FileNames fnames;
    RuntimeParameters rtp;

    // The pi 1 and 2 grids
    std::vector<long long> grid_pi1;
    std::vector<long long> grid_pi2;
    int length;

    // Results
    std::vector<long long> results1;
    std::vector<long long> results2;

    // Define the pointers
    int pointer1 = 0;
    int pointer2 = 0;

    // The outstream for this tracker
    FILE* outfile;

    // Helpers
    void _help_step_1(const long double, const long long);
    void _help_step_2(const long double, const long long);

    void _dump_outfile();

public:
    Aging(const FileNames, const RuntimeParameters);
};


class AgingConfig : public Aging
{
public:
    AgingConfig(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals);
    ~AgingConfig();
};


class AgingConfigInherentStructure : public Aging
{
public:
    AgingConfigInherentStructure(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals);
    ~AgingConfigInherentStructure();
};




class AgingBasin : public Aging
{
protected:

    std::vector<long long> vec_basin_index_1;
    std::vector<long long> vec_basin_index_2;
    std::vector<int> vec_prev_state_in_basin_1;
    std::vector<int> vec_prev_state_in_basin_2;

    double _threshold;

    RuntimeParameters rtp;
    // The basin indexes reference the last basin that a tracer was in, or the
    // basin the tracer is currently in. Whether or not a tracer is currently
    // in a basin or not can be determined by comparing the energies to the
    // thresholds in the RuntimeParameters, rtp.
    long long basin_index_1 = 0;  // First basin index
    long long basin_index_2 = 0;  // Second basin index

    // Helpers
    void _help_step_1_(const long double, const double);
    void _help_step_2_(const long double, const double);
    void _help_step(const long double, const double, const double);

    void _dump_outfile();

public:
    AgingBasin(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
};


class AgingBasinThreshold : public AgingBasin
{
public:
    AgingBasinThreshold(const FileNames, const RuntimeParameters);
    ~AgingBasinThreshold();
};

class AgingBasinAttractor : public AgingBasin
{
public:
    AgingBasinAttractor(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
    ~AgingBasinAttractor();
};

class AgingBasinThresholdInherentStructure : public AgingBasin
{
public:
    AgingBasinThresholdInherentStructure(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
    ~AgingBasinThresholdInherentStructure();
};


class AgingBasinAttractorInherentStructure : public AgingBasin
{
public:
    AgingBasinAttractorInherentStructure(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals, const Vals);
    ~AgingBasinAttractorInherentStructure();
};


#endif
