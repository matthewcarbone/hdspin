#ifndef ENERGY_MAPPING_H
#define ENERGY_MAPPING_H

#include <random>

#include "utils.h"
#include "lru.h"

class EnergyMapping
{
protected:
    utils::SimulationParameters params;

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Distributions; we'll only use one of these depending on which type
    // of dynamics we're doing
    mutable std::exponential_distribution<double> exponential_distribution;
    mutable std::normal_distribution<double> normal_distribution;

    // One must set the capacity using `set_capacity(int)`
    mutable cache::lru_cache<std::string, double> energy_map;

public:
    double sample_energy() const;
    double get_config_energy(const ap_uint<PRECISON>) const;
    void get_config_energies_array_(const ap_uint<PRECISON> *neighbors, double *neighboring_energies, const unsigned int bitLength) const;
    ap_uint<PRECISON> get_size(){return energy_map.get_size();}
    ap_uint<PRECISON> get_capacity(){return energy_map.get_capacity();}
    void _initialize_distributions();
    EnergyMapping(const utils::SimulationParameters);
    /**
     * @brief Gets the inherent structure only
     * @details [long description]
     * @return [description]
     */
    ap_uint<PRECISON> get_inherent_structure(const ap_uint<PRECISON> state) const;

};

#endif
