#ifndef SPIN_BASE_H
#define SPIN_BASE_H

#include <random>

#include "Utils/structures.h"
#include "Utils/lru.h"


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

    // Pointer to the energy mapping, only initialized if rtp.memory == -1
    double *_emap = 0;
    bool _emap_allocated = false;  // Useful flag for the destructor later

    // In the case where we're using adjusted memory dynamics
    // This is just the default constructor, it doesn't really do anything
    // One must set the capacity using `set_capacity(int)`
    mutable LRUCache lru;

public:
    double sample_energy() const;
    double get_config_energy(const long long) const;
    EnergyMapping(const RuntimeParameters);
    ~EnergyMapping();
};


class SpinSystem
{

// Accessible only from within the class or it children
protected:
    RuntimeParameters rtp;
    EnergyMapping emap;

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Pointer to the configuration in the case where we have memory
    int *_spin_config = 0;  // NULL
    bool _spin_config_allocated = false;

    // Use an index for this instead in the case when we do not have memory
    long long _memoryless_system_config;
    double _memoryless_system_energy;

    // Pointer to the inherent structure mapping
    long long *_ism = 0;
    bool _ism_allocated = false;

    // Pointer to the neighboring energies, used in the inherent structure
    // computation
    double *_neighboring_energies = 0;
    bool _neighboring_energies_allocated = true;

    // Initialize some objects for storing the previous and current values of
    // things:
    Vals prev, curr;

    // Multiplier for sampling from the waiting time and total
    // exit rate. This is 1 by default but is set to try and find the
    // equivalent Gillespie simulation for the "loop" standard dynamics.
    double _waiting_time_multiplier = 1.0;

    // Getter for the current spin representation, energies, etc.
    void _flip_spin_(const int);
    void _helper_calculate_neighboring_energies(int *, int, double *) const;
    void _calculate_neighboring_energies() const;
    long long _help_get_inherent_structure() const;
    long long _get_inherent_structure() const;

    // Updater for the previous values; this should be done at the end of
    // every recording phase
    void _init_prev();
    void _init_curr();

// Accessible outside of the class instance
public:
    SpinSystem(const RuntimeParameters, EnergyMapping);
    Vals get_prev() const {return prev;}
    Vals get_curr() const {return curr;}
    ~SpinSystem();
};

#endif
