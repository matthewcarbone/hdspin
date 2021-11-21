/* Core system utilities.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <cstring>
#include <random>

#include "Utils/structures.h"
#include "Spin/base.h"
#include "Utils/utils.h"

// Base Spin system -----------------------------------------------------------

/**
 * @brief Initializes the spin system.
 * @details Initializes the spin_config object using the Bernoulli distribution
 * generator.
 */
void BaseSpinSystem::_initialize_spin_system()
{
    // Initialize the distribution
    std::bernoulli_distribution _bernoulli_distribution;

    // Fill the vector
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        spin_config[ii] = _bernoulli_distribution(generator);
    }
}

long long BaseSpinSystem::_get_int_rep() const
{
    return binary_vector_to_int(spin_config, rtp.N_spins);
}

double BaseSpinSystem::_get_random_energy()
{
    throw std::runtime_error("_get_energy method must be overridden in inherited classes");
}

void BaseSpinSystem::init_prev_()
{
    prev.int_rep = _get_int_rep();
    prev.energy = _get_random_energy();
    prev.energy_IS = 0.0;
    prev.int_rep_IS = 0;
}

void BaseSpinSystem::init_curr_()
{
    curr.int_rep = _get_int_rep();
    curr.energy = _get_random_energy();
    curr.energy_IS = 0.0;
    curr.int_rep_IS = 0;
}

BaseSpinSystem::BaseSpinSystem(const RuntimeParameters rtp) : rtp(rtp)
{
    unsigned int seed = std::random_device{}();
    generator.seed(seed);
    spin_config = new int[rtp.N_spins];
    _initialize_spin_system();
};

BaseSpinSystem::~BaseSpinSystem(){delete[] spin_config;}




// Exponential Distribution Spin system ---------------------------------------

MemorylessExponentialSpinSystem::MemorylessExponentialSpinSystem(const RuntimeParameters rtp) : BaseSpinSystem(rtp)
{
    // Set the exponential distribution parameter
    const double p = rtp.beta_critical;
    distribution.param(std::exponential_distribution<double>::param_type(p));
}


double MemorylessExponentialSpinSystem::_get_random_energy()
{
    return -distribution(generator);
}



// Normal Distribution Spin system --------------------------------------------

MemorylessNormalSpinSystem::MemorylessNormalSpinSystem(const RuntimeParameters rtp) : BaseSpinSystem(rtp)
{
    // Set the normal distribution parameter
    const double p = sqrt(rtp.N_spins);
    distribution.param(std::normal_distribution<double>::param_type(0.0, p));
}

double MemorylessNormalSpinSystem::_get_random_energy()
{
    return distribution(generator);
}


// Auxiliary ------------------------------------------------------------------

_WithMemory::_WithMemory() {};

_WithMemory::~_WithMemory()
{
    delete[] emap;
    delete[] ism;
    delete[] neighboring_energies;
}





// Memoryless exp Spin system -------------------------------------------------


ExponentialSpinSystem::ExponentialSpinSystem(const RuntimeParameters rtp) : MemorylessExponentialSpinSystem(rtp)
{

    // Initialize the energy mapping
    emap = new double[rtp.N_configs];
    for (unsigned long long ii; ii< rtp.N_configs; ii++)
    {
        emap[ii] = _get_random_energy();
    }

    // Initialize the inherent structure mapping
    // Store every entry as -1 (to indicate that none exists yet)
    ism = new long long[rtp.N_configs];
    for (long long ii=0; ii<rtp.N_configs; ii++){ism[ii] = -1;}
}


ExponentialSpinSystem::ExponentialSpinSystem(const RuntimeParameters rtp) : MemorylessExponentialSpinSystem(rtp)
{

    // Initialize the energy mapping
    emap = new double[rtp.N_configs];
    for (unsigned long long ii; ii< rtp.N_configs; ii++)
    {
        emap[ii] = _get_random_energy();
    }

    // Initialize the inherent structure mapping
    // Store every entry as -1 (to indicate that none exists yet)
    ism = new long long[rtp.N_configs];
    for (long long ii=0; ii<rtp.N_configs; ii++){ism[ii] = -1;}
}




// void BaseSpinSystem::init_prev_()
// {
//     // Initialize the "previous" values
//     prev.int_rep = _get_int_rep();
//     prev.energy = get_energy(prev.int_rep);
//     prev.energy_IS = 0.0;
//     if (rtp.memoryless)
//     {
//         prev.energy_IS = 0.0;
//         prev.int_rep_IS = 0;
//     }
//     else
//     {
//         prev.int_rep_IS = get_inherent_structure();
//         prev.energy_IS = get_energy(prev.int_rep_IS);
//     }
// }

// void BaseSpinSystem::init_curr_()
// {
//     // Initialize the "current" values
//     curr.int_rep = get_int_rep();
//     curr.energy = get_energy(curr.int_rep);
//     if (rtp.memoryless)
//     {
//         curr.energy_IS = 0.0;
//         curr.int_rep_IS = 0;
//     }
//     else
//     {
//         curr.int_rep_IS = get_inherent_structure();
//         curr.energy_IS = get_energy(curr.int_rep_IS);
//     }
// }



void SpinSystem::_helper_calculate_neighboring_energies_(int *cfg,
    int N, double *neighboring_energies) const
{
    for (int ii=0; ii<N; ii++)
    {
        // Flip the ii'th spin
        _helper_flip_spin_(cfg, ii);

        // Collect the energy of the spin_config
        neighboring_energies[ii] =
            get_energy(binary_vector_to_int(cfg, N));

        // Flip the ii'th spin back
        _helper_flip_spin_(cfg, ii);
    }
}

void SpinSystem::_calculate_neighboring_energies()
{
    _helper_calculate_neighboring_energies_(spin_config,
        rtp.N_spins, neighboring_energies);
}


/* Finds the lowest energy configuration via greedy search from a config
 * through its neighbors. Returns the config index of this inherent
 structure. */
long long SpinSystem::_help_get_inherent_structure() const
{

    int min_el;

    // Make a copy of the config
    int config_copy[rtp.N_spins];
    memcpy(config_copy, spin_config, rtp.N_spins*sizeof(int));

    // Neighboring energies local copy
    double ne[rtp.N_spins];

    // While we are not at the inherent structure, keep going
    while (true)
    {
        // Compute the neighboring energies on the copy
        _helper_calculate_neighboring_energies_(config_copy,
            rtp.N_spins, ne);

        // If we have not reached the lowest local energy which can be reached
        // by flipping a spin
        min_el = min_element(ne, rtp.N_spins);
        if (ne[min_el] < emap[binary_vector_to_int(config_copy, rtp.N_spins)])
        {
            _helper_flip_spin_(config_copy, min_el);
        }
        else
        {
            return binary_vector_to_int(config_copy, rtp.N_spins);
        }
    }
}


void SpinSystem::flip_spin_(const int idx)
{
    _helper_flip_spin_(spin_config, idx);
}





double SpinSystem::get_current_energy() const
{
    return get_energy(get_int_rep());
}


/* Updates the inherent_structure_mapping if needed and returns the config
of the inherent structure */
long long SpinSystem::get_inherent_structure() const
{
    // Get the current configuration integer representations and energies
    const long long config_int = get_int_rep();
    
    // Update the inherent structure dictionary
    long long config_IS_int;
    if (ism[config_int] != -1)
    {
        // We've computed the inherent structure before, no need to do
        // it again
        config_IS_int = ism[config_int];
    }
    else
    {
        config_IS_int = _help_get_inherent_structure();
        ism[config_int] = config_IS_int;
    }
    return config_IS_int;
}


SpinSystem::~SpinSystem()
{
    delete[] spin_config;

    if (! rtp.memoryless)
    {
        delete[] emap;
        delete[] ism;
        delete[] neighboring_energies;
    }
}


SpinSystemWithInherentStructure::SpinSystemWithInherentStructure(
    const RuntimeParameters rtp) : SpinSystem(rtp) {};






