#include <cstring>
#include <random>

#include "Utils/structures.h"
#include "Spin/base.h"
#include "Utils/utils.h"

// Base Spin system -----------------------------------------------------------

SpinSystem::SpinSystem(const RuntimeParameters rtp) : rtp(rtp)
{
    unsigned int seed = std::random_device{}();
    generator.seed(seed);

    // Initialize the distributions themselves
    if (rtp.landscape == 0)
    {
        const double p = rtp.beta_critical;
        exponential_distribution.param(std::exponential_distribution<double>::param_type(p));
    }
    else
    {
        const double p = sqrt(rtp.N_spins);
        normal_distribution.param(std::normal_distribution<double>::param_type(0.0, p));
    }

    // Allocate memory for the spin configuration
    spin_config = new int[rtp.N_spins];

    // Initialize the distribution
    std::bernoulli_distribution _bernoulli_distribution;

    // Fill the spin_config vector
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        spin_config[ii] = _bernoulli_distribution(generator);
    }

    // If we do not have a memoryless system, initialize the appropriate
    // memory
    if (rtp.memoryless == 0)
    {
        // Initialize the energy mapping
        emap = new double[rtp.N_configs];
        for (unsigned long long ii=0; ii<rtp.N_configs; ii++)
        {
            emap[ii] = _get_random_energy();
        }

        // Initialize the inherent structure mapping
        // Store every entry as -1 (to indicate that none exists yet)
        ism = new long long[rtp.N_configs];
        for (long long ii=0; ii<rtp.N_configs; ii++){ism[ii] = -1;}
        spin_config_energy = -999999.0;
    }
    else
    {
        spin_config_energy = _get_random_energy();
    }

    if (rtp.dynamics_flag > 1){_waiting_time_multiplier = 1.0 / rtp.N_spins;}
};

/**
 * @brief Flips a spin
 * @details Flips the spin at site idx
 *
 * @param idx Location to flip the spin
 */
void SpinSystem::_flip_spin_(const int idx)
{
    _helper_flip_spin_(spin_config, idx);
}

/**
 * @details Gets the integer representation of the spin_config vector
 */
long long SpinSystem::_get_int_rep() const
{
    return binary_vector_to_int(spin_config, rtp.N_spins);
}

/**
 * @details Depending on the landscape, returns a random sample from the
 * appropriate distribution
 */
double SpinSystem::_get_random_energy() const
{
    if (rtp.landscape == 0){return -exponential_distribution(generator);}
    return normal_distribution(generator);
}

/**
 * @details Gets the energy corresponding to the integer representation
 * provided
 */
double SpinSystem::_get_energy(const long long int_rep) const
{
    if (rtp.memoryless == 1){return spin_config_energy;}
    return emap[int_rep];
}

double SpinSystem::_get_current_energy() const
{
    return _get_energy(binary_vector_to_int(spin_config, rtp.N_spins));
}

void SpinSystem::_helper_calculate_neighboring_energies(int *cfg,
    int N, double *ne) const
{
    for (int ii=0; ii<N; ii++)
    {
        if (rtp.memoryless == 0)
        {
            // Flip the ii'th spin
            _helper_flip_spin_(cfg, ii);

            // Collect the energy of the spin_config
            ne[ii] = _get_energy(binary_vector_to_int(cfg, N));

            // Flip the ii'th spin back
            _helper_flip_spin_(cfg, ii);
        }
        else
        {
            ne[ii] = _get_random_energy();
        }
    }
}

void SpinSystem::_calculate_neighboring_energies() const
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
        _helper_calculate_neighboring_energies(config_copy, rtp.N_spins, ne);

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

    throw std::runtime_error("Unknown error in _help_get_inherent_structure");
}


/* Updates the inherent_structure_mapping if needed and returns the config
of the inherent structure */
long long SpinSystem::_get_inherent_structure() const
{
    if (rtp.memoryless == 1)
    {
        throw std::runtime_error("Inherent structure not defined for memoryless spin systems");
    }

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


void SpinSystem::_init_prev()
{
    prev.int_rep = _get_int_rep();
    prev.energy = _get_energy(prev.int_rep);

    if (rtp.memoryless == 0) // Using memory
    {
        prev.int_rep_IS = _get_inherent_structure();
        prev.energy_IS = _get_energy(prev.int_rep_IS);
    }
}

void SpinSystem::_init_curr()
{
    curr.int_rep = _get_int_rep();
    curr.energy = _get_energy(curr.int_rep);
    if (rtp.memoryless == 0) // Using memory
    {
        curr.int_rep_IS = _get_inherent_structure();
        curr.energy_IS = _get_energy(curr.int_rep_IS);
    }
}


SpinSystem::~SpinSystem()
{
    delete[] spin_config;
    if (rtp.memoryless == 0) // Using memory
    {
        delete[] emap;
        delete[] ism;
    }
}



