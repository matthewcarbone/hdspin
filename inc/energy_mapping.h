#ifndef ENERGYMAPPING_H
#define ENERGYMAPPING_H

#include <random>

#include "structures.h"
#include "lru.h"

class EnergyMapping
{
protected:
    RuntimeParameters rtp;

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Distributions; we'll only use one of these depending on which type
    // of dynamics we're doing
    mutable std::exponential_distribution<double> exponential_distribution;
    mutable std::normal_distribution<double> normal_distribution;

    // One must set the capacity using `set_capacity(int)`
    mutable LRUCache energy_map;

public:
    double sample_energy() const;
    double get_config_energy(const long long) const;
    EnergyMapping(const RuntimeParameters);
    ~EnergyMapping();
};


#endif
