#include <cstring>
#include <random>

#include "Utils/structures.h"
#include "Spin/base.h"
#include "Utils/utils.h"

// Base Spin system -----------------------------------------------------------

MemorylessSpinSystem::MemorylessSpinSystem(const RuntimeParameters rtp) : rtp(rtp)
{
    unsigned int seed = std::random_device{}();
    generator.seed(seed);
    spin_config = new int[rtp.N_spins];
    _initialize_spin_system();
};

void MemorylessSpinSystem::_flip_spin(const int idx)
{
    _helper_flip_spin_(spin_config, idx);
}


/**
 * @brief Initializes the spin system.
 * @details Initializes the spin_config object using the Bernoulli distribution
 * generator.
 */
void MemorylessSpinSystem::_initialize_spin_system()
{
    // Initialize the distribution
    std::bernoulli_distribution _bernoulli_distribution;

    // Fill the vector
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        spin_config[ii] = _bernoulli_distribution(generator);
    }
}

long long MemorylessSpinSystem::_get_int_rep() const
{
    return binary_vector_to_int(spin_config, rtp.N_spins);
}

double MemorylessSpinSystem::_get_random_energy() const
{
    if (rtp.landscape == 0){return -exponential_distribution(generator);}
    else {return normal_distribution(generator);}
}

double MemorylessSpinSystem::_get_energy(const long long int_rep) const
{
    return _get_random_energy();
}

long long MemorylessSpinSystem::get_inherent_structure() const
{
    throw std::runtime_error("Inherent structure not defined for memoryless spin systems");
}

void MemorylessSpinSystem::init_prev_()
{
    prev.int_rep = _get_int_rep();
    prev.energy = _get_random_energy();
}

void MemorylessSpinSystem::init_curr_()
{
    curr.int_rep = _get_int_rep();
    curr.energy = _get_random_energy();
}

MemorylessSpinSystem::~MemorylessSpinSystem(){delete[] spin_config;}



// Auxiliary ------------------------------------------------------------------

SpinSystem::SpinSystem(const RuntimeParameters rtp) : MemorylessSpinSystem(rtp)
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
};

double SpinSystem::_get_energy(const long long int_rep) const
{
    return emap[int_rep];
}

void SpinSystem::init_prev_()
{
    prev.int_rep = _get_int_rep();
    prev.energy = _get_energy(prev.int_rep);
    prev.int_rep_IS = get_inherent_structure();
    prev.energy_IS = _get_energy(prev.int_rep_IS);
}

void SpinSystem::init_curr_()
{
    curr.int_rep = _get_int_rep();
    curr.energy = _get_energy(curr.int_rep);
    curr.int_rep_IS = get_inherent_structure();
    curr.energy_IS = _get_energy(curr.int_rep_IS);
}


void SpinSystem::_helper_calculate_neighboring_energies(int *cfg,
    int N, double *neighboring_energies) const
{
    for (int ii=0; ii<N; ii++)
    {
        // Flip the ii'th spin
        _helper_flip_spin_(cfg, ii);

        // Collect the energy of the spin_config
        neighboring_energies[ii] = _get_energy(binary_vector_to_int(cfg, N));

        // Flip the ii'th spin back
        _helper_flip_spin_(cfg, ii);
    }
}

void SpinSystem::_calculate_neighboring_energies()
{
    _helper_calculate_neighboring_energies(spin_config,
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
        _helper_calculate_neighboring_energies(config_copy,
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

/* Updates the inherent_structure_mapping if needed and returns the config
of the inherent structure */
long long SpinSystem::get_inherent_structure() const
{
    // Get the current configuration integer representations and energies
    const long long config_int = _get_int_rep();

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
    delete[] emap;
    delete[] ism;
    delete[] neighboring_energies;
}
