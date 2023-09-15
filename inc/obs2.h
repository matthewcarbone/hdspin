#ifndef OBS_PSI_H
#define OBS_PSI_H

#include <vector>

#include "spin.h"
#include "utils.h"


class PsiConfig
{
protected:

    const parameters::FileNames fnames;
    const parameters::SimulationParameters params;
    FILE *outfile;
    const SpinSystem* spin_system_ptr;
    
    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    long long _max_counter;

    std::vector<long long> _counter;

    // Keep track internally of the waiting time for both the standard
    // trajectory and the inherent structure
    long double _waiting_time = 0.0;

    int _out_of_counter = 0;

public:

    PsiConfig(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double current_waiting_time);
    ~PsiConfig();

};



#endif
