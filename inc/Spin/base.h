#ifndef SPIN_BASE_H
#define SPIN_BASE_H

#include <random>

#include "Utils/structures.h"


class SpinSystem
{

// Accessible only from within the class or it children
protected:
    RuntimeParameters rtp;

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Distributions; we'll only use one of these depending on which type
    // of dynamics we're doing
    mutable std::exponential_distribution<double> exponential_distribution;
    mutable std::normal_distribution<double> normal_distribution;

    // Pointer to the configuration
    int *spin_config = 0;  // NULL

    // Pointer to the energy mapping
    double *emap = 0;

    // Pointer to the inherent structure mapping
    long long *ism = 0;

    // Pointer to the neighboring energies, used in the inherent structure
    // computation
    double *neighboring_energies = 0;

    // Initialize some objects for storing the previous and current values of
    // things:
    Vals prev, curr;

    // Number of accepted steps (non-rejections)
    long long n_accept = 0;

    // Multiplier for sampling from the waiting time and total
    // exit rate. This is 1 by default but is set to try and find the
    // equivalent Gillespie simulation for the "loop" standard dynamics.
    double _waiting_time_multiplier = 1.0;

    // Getter for the current spin representation, energies, etc.
    void _flip_spin_(const int);
    long long _get_int_rep() const;
    double _get_random_energy() const;
    double _get_energy(const long long) const;
    double _get_current_energy() const;
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
    SpinSystem(const RuntimeParameters);
    Vals get_prev() const {return prev;}
    Vals get_curr() const {return curr;}
    int get_n_accept() const {return n_accept;}
    ~SpinSystem();
};

#endif
