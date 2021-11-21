#ifndef SIM_H
#define SIM_H

#include <random>

#include "Utils/structures.h"


class Simulation
{
protected:
    // Global parameters for the simulation
    FileNames fnames;
    RuntimeParameters rtp;

    // Simulation parameters
    long double simulation_clock = 0.0;

public:
    Simulation(const FileNames, const RuntimeParameters);

};


class GillespieSimulation : public Simulation
{
public:
    GillespieSimulation(const FileNames, const RuntimeParameters);
    void execute();
};

class GillespieSimulationMemoryless : public Simulation
{
public:
    GillespieSimulationMemoryless(const FileNames, const RuntimeParameters);
    void execute();
};

class StandardSimulation : public Simulation
{
public:
    StandardSimulation(const FileNames, const RuntimeParameters);
    void execute();
};

class StandardSimulationMemoryless : public Simulation
{
public:
    StandardSimulationMemoryless(const FileNames, const RuntimeParameters);
    void execute();
};

#endif
