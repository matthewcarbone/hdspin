#ifndef SPIN_BASE_H
#define SPIN_BASE_H

#include <random>

#include "Utils/structures.h"


class SpinSystem
{
protected:
    RuntimeParameters rtp;

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    std::mt19937 generator;

    // Pointer to the configuration
    int *spin_config = 0;  // NULL

    // Pointer to the energy mapping
    double *emap = 0;

    // Pointer to the inherent structure mapping
    long long *ism = 0;

    // Pointer to the neighboring energies, used in the inherent structure
    // computation
    double *neighboring_energies = 0;

    // Number of accepted steps (non-rejections)
    long long n_accept = 0;

    // Initialize some objects for storing the previous and current values of
    // things:
    Vals prev, curr;

    // Initializes the spin system by clearing the current configuration and
    // filling it with new, randomly selected up/down (1/0) binary values.
    // This is called once at instantiation.
    void _initialize_spin_system();
    void _initialize_energy_mapping_();
    void _initialize_inherent_structure_mapping_();

    // Updater for the previous values; this should be done at the end of
    // every recording phase
    void init_prev_();
    void init_curr_();

    void _calculate_neighboring_energies();
    long long _help_get_inherent_structure() const;

public:
    SpinSystem(const RuntimeParameters);

    // Flips the spin at the specified location.
    void flip_spin_(const int);

    // Getters
    long long get_int_rep() const;
    double get_current_energy() const;
    long long get_inherent_structure() const;
    std::vector<int> get_spin_config() const;
    Vals get_prev() const {return prev;}
    Vals get_curr() const {return curr;}
    long long get_n_accept() const {return n_accept;}

    // Release memory in the destructor
    ~SpinSystem();
};


#endif
