/* Core system utilities.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <random>

#include "Utils/structures.h"
#include "Spin/base.h"
#include "Utils/utils.h"


void SpinSystem::_initialize_spin_system()
{
    // Initialize the distribution
    std::bernoulli_distribution distribution;
    
    // Fill the vector
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        spin_config[ii] = distribution(generator);
    }
}

void SpinSystem::_initialize_energy_mapping_()
{
    // Initialize prior with
    // double *energy_arr = new double[rtp.N_configs];

    if (rtp.landscape == 0)
    {
        // C++ exponential random number generator has input of the form
        // lambda * e^{-lambda * x} for x > 0. For us, the scale parameter is
        // bc.
        std::exponential_distribution<double> distribution(rtp.beta_critical);

        // Initialize the exponential distribution
        for (unsigned long long ii=0; ii<rtp.N_configs; ii++)
        {
            // Note for the EREM model the energies are all negative, hence why
            // we take the negative sign here.
            emap[ii] = -distribution(generator);
        }
    }

    else if (rtp.landscape == 1)
    {
        std::normal_distribution<double> distribution(0.0, sqrt(rtp.N_spins));

        for (unsigned long long ii=0; ii<rtp.N_configs; ii++)
        {
            emap[ii] = distribution(generator);
        }
    }

    else{throw std::runtime_error("Unknown landscape");}
}


void SpinSystem::_initialize_inherent_structure_mapping_()
{
    // Store every entry as -1 (to indicate that none exists yet)
    for (long long ii=0; ii<rtp.N_configs; ii++){ism[ii] = -1;}
}


SpinSystem::SpinSystem(const RuntimeParameters rtp) : rtp(rtp)
{
    unsigned int seed = std::random_device{}();
    generator.seed(seed);

    spin_config = new int[rtp.N_spins];
    _initialize_spin_system();
    emap = new double[rtp.N_configs];
    _initialize_energy_mapping_();
    ism = new long long[rtp.N_configs];
    _initialize_inherent_structure_mapping_();
}


void SpinSystem::init_prev_()
{
    // Initialize the "previous" values
    prev.int_rep = get_int_rep();
    prev.int_rep_IS = get_inherent_structure();
    prev.energy = emap[prev.int_rep];
    prev.energy_IS = emap[prev.int_rep_IS];
}

void SpinSystem::init_curr_()
{
    // Initialize the "current" values
    curr.int_rep = get_int_rep();
    curr.int_rep_IS = get_inherent_structure();
    curr.energy = emap[curr.int_rep];
    curr.energy_IS = emap[curr.int_rep_IS];
}

void SpinSystem::_calculate_neighboring_energies()
{
    _helper_calculate_neighboring_energies_(spin_config, emap,
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
        _helper_calculate_neighboring_energies_(config_copy, emap,
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


long long SpinSystem::get_int_rep() const
{
    return binary_vector_to_int(spin_config, rtp.N_spins);
}


double SpinSystem::get_current_energy() const
{
    return emap[get_int_rep()];
}


std::vector<int> SpinSystem::get_spin_config() const
{
    std::vector<int> v;
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        v.push_back(spin_config[ii]);
    }

    return v;
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
    delete[] emap;
    delete[] ism;
    delete[] neighboring_energies;
}
