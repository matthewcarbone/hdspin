#include <cstring>
#include <random>

#include "Utils/structures.h"
#include "Spin/base.h"
#include "Utils/utils.h"
#include "Utils/lru.h"


// Energy mapping -------------------------------------------------------------

double EnergyMapping::sample_energy() const
{
    if (rtp.landscape == "EREM"){return -exponential_distribution(generator);}
    else{return normal_distribution(generator);}
}


double EnergyMapping::get_config_energy(const long long int_rep) const
{

    // First case, we're saving all results cached in an array
    if (rtp.memory == -1){return _emap[int_rep];}

    // Second case, we're using adjusted memory dynamics, only caching some
    // finite number of energies
    else if (rtp.memory > 0)
    {
        // If our key exists in the LRU cache, simply return the value
        if (lru.key_exists(int_rep)){return lru.get_fast(int_rep);}

        // Otherwise, we sample a new value
        else
        {
            const double sampled = sample_energy();
            lru.put(int_rep, sampled);
            return sampled;
        }
    }

    else{throw std::runtime_error(
        "Do not call get_config_energy in memoryless system");}
}


EnergyMapping::EnergyMapping(const RuntimeParameters rtp) : rtp(rtp)
{
    unsigned int seed = std::random_device{}();
    generator.seed(seed);

    // Initialize the distributions themselves
    if (rtp.landscape == "EREM")
    {
        const double p = rtp.beta_critical;
        exponential_distribution.param(std::exponential_distribution<double>::param_type(p));
    }
    else
    {
        const double p = sqrt(rtp.N_spins);
        normal_distribution.param(
            std::normal_distribution<double>::param_type(0.0, p));
    }

    // Allocate the energy storage mediums
    if (rtp.memory == -1)  // Use emap directly
    {
        // Initialize the energy mapping
        _emap = new double[rtp.N_configs];
        _emap_allocated = true;
        for (long long ii=0; ii<rtp.N_configs; ii++)
        {
            _emap[ii] = sample_energy();
        }
    }
    else if (rtp.memory > 0)  // LRU queue!
    {
        lru.set_capacity(rtp.memory);
    }
};

EnergyMapping::~EnergyMapping()
{
    if (_emap_allocated){delete[] _emap;}
}


// Base Spin system -----------------------------------------------------------

SpinSystem::SpinSystem(const RuntimeParameters rtp, EnergyMapping emap)
    : rtp(rtp), emap(emap)
{
    unsigned int seed = std::random_device{}();
    generator.seed(seed);

    // In the case where we have any sort of memory, fill the spin config
    // object
    if (rtp.memory != 0)
    {
        // Allocate memory for the spin configuration
        _spin_config = new int[rtp.N_spins];
        _spin_config_allocated = true;

        // Initialize the distribution
        std::bernoulli_distribution _bernoulli_distribution;

        for (int ii=0; ii<rtp.N_spins; ii++)
        {
            _spin_config[ii] = _bernoulli_distribution(generator);
        }

        // Allocate the inherent structure mapping
        _ism = new long long[rtp.N_configs];
        _ism_allocated = true;
        for (long long ii=0; ii<rtp.N_configs; ii++){_ism[ii] = -1;}
    }
    else
    {
        _memoryless_system_config = 1;
        _memoryless_system_energy = emap.sample_energy();
    }

    if (rtp.divN == 1){_waiting_time_multiplier = 1.0 / rtp.N_spins;}
};

/**
 * @brief Flips a spin
 * @details Flips the spin at site idx
 *
 * @param idx Location to flip the spin
 */
void SpinSystem::_flip_spin_(const int idx)
{
    if (rtp.memory == 0)
    {
        throw std::runtime_error("Cannot flip spin in memoryless system");
    }
    _helper_flip_spin_(_spin_config, idx);
}

void SpinSystem::_helper_calculate_neighboring_energies(int *cfg,
    int N, double *ne) const
{
    for (int ii=0; ii<N; ii++)
    {
        if (rtp.memory != 0)
        {
            // Flip the ii'th spin
            _helper_flip_spin_(cfg, ii);

            // Collect the energy of the spin_config
            ne[ii] = emap.get_config_energy(binary_vector_to_int(cfg, N));

            // Flip the ii'th spin back
            _helper_flip_spin_(cfg, ii);
        }
        else{ne[ii] = emap.sample_energy();}
    }
}

void SpinSystem::_calculate_neighboring_energies() const
{
    _helper_calculate_neighboring_energies(_spin_config,
        rtp.N_spins, _neighboring_energies);
}

/* Finds the lowest energy configuration via greedy search from a config
 * through its neighbors. Returns the config index of this inherent
 structure. */
long long SpinSystem::_help_get_inherent_structure() const
{

    int min_el;
    double tmp_energy;

    // Make a copy of the config
    int config_copy[rtp.N_spins];
    memcpy(config_copy, _spin_config, rtp.N_spins*sizeof(int));

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
        tmp_energy = emap.get_config_energy(
            binary_vector_to_int(config_copy, rtp.N_spins));

        if (ne[min_el] < tmp_energy)
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
    if (rtp.memory == 0)
    {
        throw std::runtime_error(
            "Inherent structure not defined for memoryless spin systems");
    }

    // Get the current configuration integer representations and energies
    const long long config_int
        = binary_vector_to_int(_spin_config, rtp.N_spins);

    // Update the inherent structure dictionary
    long long config_IS_int;
    if (_ism[config_int] != -1)
    {
        // We've computed the inherent structure before, no need to do
        // it again
        config_IS_int = _ism[config_int];
    }
    else
    {
        config_IS_int = _help_get_inherent_structure();
        _ism[config_int] = config_IS_int;
    }
    return config_IS_int;
}


void SpinSystem::_init_prev()
{
    prev.int_rep = binary_vector_to_int(_spin_config, rtp.N_spins);
    prev.energy = emap.get_config_energy(prev.int_rep);
    if (rtp.memory != 0)
    {
        prev.int_rep_IS = _get_inherent_structure();
        prev.energy_IS = emap.get_config_energy(prev.int_rep_IS);
    }
}

void SpinSystem::_init_curr()
{
    curr.int_rep = binary_vector_to_int(_spin_config, rtp.N_spins);
    curr.energy = emap.get_config_energy(curr.int_rep);
    if (rtp.memory != 0)
    {
        curr.int_rep_IS = _get_inherent_structure();
        curr.energy_IS = emap.get_config_energy(curr.int_rep_IS);
    }
}


SpinSystem::~SpinSystem()
{
    if (_spin_config_allocated){delete[] _spin_config;}
    if (_ism_allocated){delete[] _ism;}
    if (_neighboring_energies_allocated){delete[] _neighboring_energies;}
}



