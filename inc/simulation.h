#ifndef SIM_H
#define SIM_H

#include <random>

#include "utils.h"


class Simulation
{
protected:
    // Global parameters for the simulation
    parameters::FileNames fnames;
    parameters::SimulationParameters params;

public:
    Simulation(const parameters::FileNames, const parameters::SimulationParameters);
    void execute();
};

#endif
