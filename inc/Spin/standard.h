#ifndef SPIN_STANDARD_H
#define SPIN_STANDARD_H

#include <random>

#include "Utils/structures.h"
#include "Spin/base.h"


class StandardSpinSystem : public SpinSystem
{
private:
    double time_in_config = 1.0;  // Minimum of 1 timestep in configuration
    double time_in_config_IS = 1.0;

public:
    StandardSpinSystem(const RuntimeParameters);

    // Step executes a possible alteration in the state, but not always. Thus,
    // the standard step actually returns whether or not the new state was
    // accepted: if there was a rejection, return false, else, if the proposed
    // state was accepted, return true.
    bool step_();
};


#endif
