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

    // Simulation parameters
    long double simulation_clock = 0.0;

public:
    Simulation(const parameters::FileNames, const parameters::SimulationParameters);
};


class GillespieSimulation : public Simulation
{
public:
    GillespieSimulation(const parameters::FileNames, const parameters::SimulationParameters);
    void execute();
};

class StandardSimulation : public Simulation
{
public:
    StandardSimulation(const parameters::FileNames, const parameters::SimulationParameters);
    void execute();
};

#endif
