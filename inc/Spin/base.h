#ifndef SPIN_BASE_H
#define SPIN_BASE_H

#include <random>

#include "Utils/structures.h"


class BaseSpinSystem
{
protected:
    RuntimeParameters rtp;

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Pointer to the configuration
    int *spin_config = 0;  // NULL

    // Initialize some objects for storing the previous and current values of
    // things:
    Vals prev, curr;

    // Number of accepted steps (non-rejections)
    long long n_accept = 0;

    // Getter for the current spin representation
    long long _get_int_rep() const;

    // Updater for the previous values; this should be done at the end of
    // every recording phase
    void init_prev_();
    void init_curr_();

    virtual double _get_random_energy();

    // Initializes spin_config
    void _initialize_spin_system();

public:
    BaseSpinSystem(const RuntimeParameters);
    ~BaseSpinSystem();
};


class MemorylessExponentialSpinSystem : public BaseSpinSystem
{
protected:
    std::exponential_distribution<double> distribution;
public:
    MemorylessExponentialSpinSystem(const RuntimeParameters);
    double _get_random_energy();
};


class MemorylessNormalSpinSystem : public BaseSpinSystem
{
protected:
    std::normal_distribution<double> distribution;
public:
    MemorylessNormalSpinSystem(const RuntimeParameters);
    double _get_random_energy();
};


class _WithMemory
{
protected:
    // Pointer to the energy mapping
    double *emap = 0;

    // Pointer to the inherent structure mapping
    long long *ism = 0;

    // Pointer to the neighboring energies, used in the inherent structure
    // computation
    double *neighboring_energies = 0;

public:
    _WithMemory();
    ~_WithMemory();
};



// Note that the inherent structure is not defined for systems without memory
class ExponentialSpinSystem : public MemorylessExponentialSpinSystem, public _WithMemory
{
protected:

    // Initializers for the above
    void _initialize_inherent_structure_mapping_();

public:
    ExponentialSpinSystem(const RuntimeParameters);
};



class MemorylessSpinSystem
{
protected:


    // Pointer to the energy mapping
    double *emap = 0;

    // Number of accepted steps (non-rejections)
    long long n_accept = 0;

    // Initializes the spin system by clearing the current configuration and
    // filling it with new, randomly selected up/down (1/0) binary values.
    // This is called once at instantiation.
    void _initialize_spin_system();
    void _initialize_energy_mapping_();



public:
    MemorylessSpinSystem(const RuntimeParameters);

    // Flips the spin at the specified location.
    void flip_spin_(const int);

    // Getters
    long long get_int_rep() const;
    double get_current_energy() const;
    double get_energy(const long long) const;
    long long get_inherent_structure() const;
    std::vector<int> get_spin_config() const;
    Vals get_prev() const {return prev;}
    Vals get_curr() const {return curr;}
    long long get_n_accept() const {return n_accept;}

    // Release memory in the destructor
    ~MemorylessSpinSystem();
};


class SpinSystemWithInherentStructure : public MemorylessSpinSystem
{
protected:

    // Pointer to the inherent structure mapping
    long long *ism = 0;

    // Pointer to the neighboring energies, used in the inherent structure
    // computation
    double *neighboring_energies = 0;

    void _initialize_inherent_structure_mapping_();
    void _calculate_neighboring_energies();
    long long _help_get_inherent_structure() const;

    /* Computes the neighboring energies of a configuration. Takes as input the
    configuration, energy array and length of the configuration, and fills in the
    fourth argument, the neighboring_energies, with the energies of the neighbors
    acquired by flipping that respective spin. */
    void _helper_calculate_neighboring_energies_(int *, int, double *) const;

public:
    SpinSystemWithInherentStructure(const RuntimeParameters);
    ~SpinSystemWithInherentStructure();
};


#endif
