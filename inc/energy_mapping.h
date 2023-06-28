#ifndef ENERGYMAPPING_H
#define ENERGYMAPPING_H

#include <random>

#include "utils.h"
#include "lru.h"

class EnergyMapping
{
protected:
    parameters::SimulationParameters parameters;

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
    double get_config_energy(const ap_uint<PRECISON>) const;
    void get_config_energies_array(const ap_uint<PRECISON> *neighbors, double *neighboring_energies, const int bitLength);
    void _initialize_distributions();
    EnergyMapping(const parameters::SimulationParameters);
};

#endif
