#ifndef OBS_AGE_H
#define OBS_AGE_H

#include <vector>

#include "Utils/structures.h"

class AgingConfig
{
private:

    // The pi 1 and 2 grids
    std::vector<long long> grid_pi1;
    std::vector<long long> grid_pi2;
    int length_pi1, length_pi2;
    long long max_time_pi1, max_time_pi2;

    // Define the pointers
    int pointer1 = 0;
    int pointer2 = 0;

    // The outstream for this tracker
    FILE *outfile_pi1;
    FILE *outfile_pi2;

    // Helpers
    void _help_step_1_(const long double, const long long, const Vals);
    void _help_step_2_(const long double, const long long, const Vals);

public:

    AgingConfig(const FileNames);
    void step_(const long double, const long long, const Vals);
    ~AgingConfig();
};


#endif
