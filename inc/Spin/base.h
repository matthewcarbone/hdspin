#ifndef SPIN_BASE_H
#define SPIN_BASE_H

#include <random>

#include "Utils/structures.h"


class MemorylessSpinSystem
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

    // Pointer to the configuration
    int *spin_config = 0;  // NULL

    // Initialize some objects for storing the previous and current values of
    // things:
    Vals prev, curr;

    // Number of accepted steps (non-rejections)
    long long n_accept = 0;

    // Getter for the current spin representation, energies, etc.
    void _flip_spin(const int);
    long long _get_int_rep() const;
    double _get_random_energy() const;
    virtual double _get_energy(const long long) const;
    virtual long long get_inherent_structure() const;

    // Updater for the previous values; this should be done at the end of
    // every recording phase
    void init_prev_();
    void init_curr_();

    // Initializes spin_config
    void _initialize_spin_system();

public:
    MemorylessSpinSystem(const RuntimeParameters);
    Vals get_prev() const {return prev;}
    Vals get_curr() const {return curr;}
    ~MemorylessSpinSystem();
};



class SpinSystem : public MemorylessSpinSystem
{
protected:
    // Pointer to the energy mapping
    double *emap = 0;

    // Pointer to the inherent structure mapping
    long long *ism = 0;

    // Pointer to the neighboring energies, used in the inherent structure
    // computation
    double *neighboring_energies = 0;

    // Getter for the current spin representation, energies, etc.
    double _get_energy(const long long) const;

    // Updater for the previous values; this should be done at the end of
    // every recording phase
    void init_prev_();
    void init_curr_();

    /* Computes the neighboring energies of a configuration. Takes as input the
    configuration, energy array and length of the configuration, and fills in the
    fourth argument, the neighboring_energies, with the energies of the neighbors
    acquired by flipping that respective spin. */
    void _helper_calculate_neighboring_energies(int *, int, double *) const;
    void _calculate_neighboring_energies();
    long long _help_get_inherent_structure() const;
    long long get_inherent_structure() const;

public:
    SpinSystem(const RuntimeParameters);
    ~SpinSystem();
};

#endif
