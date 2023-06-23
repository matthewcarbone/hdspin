#ifndef CONFIGS_H
#define CONFIGS_H

#include <unordered_set>
#include <vector>

#include "Utils/structures.h"


class Configs
{
protected:
    FileNames fnames;
    RuntimeParameters rtp;
    FILE* outfile;
    std::unordered_set<int> _unique_configs;
    // The maximum time on the grid
    long long max_time;

    // The grid (in time) itself
    std::vector<long long> grid;

    // The length of the grid
    int length;

    // The pointer to the last-updated point on the grid
    int pointer = 0;

public:
    Configs(const FileNames, const RuntimeParameters);
    void step(const long double, const Vals);
    ~Configs();
};

#endif
