#ifndef ENERGY_MAPPING_H
#define ENERGY_MAPPING_H

#include <random>

#include "utils.h"
#include "lru.h"

/**
 * @brief Provides the apparatus for a mapping between configurations of spins
 * and their energies
 * 
 * @details The EnergyMapping class uses a least-recently used queue for storing
 * energy values, with string representations of AP integers as the keys (the
 * AP integer type isn't hashable). It also contains the distributional information
 * for sampling energies from either the exponential distribution or the normal
 * distribution. These are parameterized by beta and beta critical, which the
 * user provides (to be clear, beta is provided, beta critical is a function of the
 * system).
 * 
 * The class comes with a few helper methods. One can sample a random energy via
 * sample_energy, get the energy of a configuration by providing its AP integer
 * representation to get_config_energy, or get all of the neighbor spin state
 * energies via get_config_energies_array_.
 * 
 */
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

    // Self-explanatory: initializes one of the above two distributions
    // depending on the simulation parameters.
    void _initialize_distributions();

    // One must set the capacity using `set_capacity(int)`, since this is the
    // zero-parameter constructor (the only type allowed in this case)
    mutable cache::lru_cache<std::string, double> energy_map;

public:

    /**
     * @brief Samples energy randomly
     * @details Samples energy randomly using the distribution specified by the
     * input parameters in the constructor.
     * @return The sampled energy
     */
    double sample_energy() const;

    /**
     * @brief Gets the energy of a provided configuration
     * @details Given an AP integer, gets the energy by converting the provided
     * integer into its string representation and querying the energy_map. If no
     * key exists, it is sampled.
     * 
     * @param state The AP integer representation of the state
     * @return The energy of that configuration
     */
    double get_config_energy(const ap_uint<PRECISON> state) const;

    /**
     * @brief Gets every neighboring configuration's energy given a config
     * @details Given some state [s1, s2, s3, ..., sN], where si are the {0, 1}
     * binary spins, this function fills a provided array (neighboring_energies)
     * with the energy of the configuration reached by flipping the ith bit. In
     * other words, neighboring_energies[ii] is the energy of the state after bit
     * ii has been flipped (counting left to right).
     * 
     * @param neighbors The neighboring sites
     * @param neighboring_energies The neighboring energies to be filled
     * @param int The total number of spins aka the bit length
     */
    void get_config_energies_array_(const ap_uint<PRECISON> *neighbors, double *neighboring_energies, const unsigned int bitLength) const;

    /**
     * @brief Returns the current size of the LRU cache
     */
    ap_uint<PRECISON> get_size() const {return energy_map.get_size();}

    /**
     * @brief Gets the maximum capacity of the LRU cache.
     */
    ap_uint<PRECISON> get_capacity() const {return energy_map.get_capacity();}

    /**
     * @brief Gets the inherent structure of the state
     * @details The inherent structure is defined as the structure which can be
     * greedily reached by flipping one spin at a time, such that flipping the
     * ith spin takes one to a lower energy. Once this is no longer possible, that is the inherent structure.
     * @return The AP representation of the inherent structure
     */
    ap_uint<PRECISON> get_inherent_structure(const ap_uint<PRECISON> state) const;
    
    /**
     * @brief Constructor for EnergyMapping
     * 
     * @param The parameters of the simulation
     */
    EnergyMapping(const utils::SimulationParameters);
    

};

#endif
