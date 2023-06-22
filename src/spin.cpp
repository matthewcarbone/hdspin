#include <cstring>
#include <random>

#include "utils.h"
#include "spin.h"


void SpinSystem::_initialize_state_()
{
    if (params.use_manual_seed == true)
    {
        generator.seed(params.seed);
    }
    else
    {
        const unsigned int seed = std::random_device{}();  
        generator.seed(seed);
    }

    int* spin_config = 0;
    spin_config = new int [params.N_spins];

    std::bernoulli_distribution _bernoulli_distribution;

    for (int ii=0; ii<params.N_spins; ii++)
    {
        spin_config[ii] = _bernoulli_distribution(generator);
    }

    state::arbitrary_precision_integer_from_int_array_(spin_config, params.N_spins, current_state);

    delete[] spin_config;
}

void SpinSystem::_fill_neighbors_()
{
    state::get_neighbors_(neighbors, current_state, params.N_spins);
}

void SpinSystem::_fill_neighboring_energies_()
{
    for (int ii=0; ii<params.N_spins; ii++)
    {
        neighboring_energies[ii] = emap_ptr->get_config_energy(neighbors[ii]);
    }
}


SpinSystem::SpinSystem(const parameters::SimulationParameters params,
    EnergyMapping& emap) : params(params)
{
    emap_ptr = &emap;

    _initialize_state_();

    // std::uniform_int_distribution<> int_dist(0, params.N_spins - 1);

    // if (params.divN == 1){_waiting_time_multiplier = 1.0 / params.N_spins;}

    // Initialize an array for containing the neighbors themselves
    neighbors = new ap_uint<PRECISON>[params.N_spins];

    // Initialize an array for containing the neighboring energies of the
    // state
    neighboring_energies = new double[params.N_spins];

    // _spin_config if in the memoryless calculation will just be a null
    // pointer; nothing is done in the helper to the config if we're doing
    // memoryless, so this should work fine
    // _helper_fill_neighboring_energies(_spin_config, params.N_spins, _neighboring_energies);

    // Initialize the previous energy to something random, from one of
    // the neighbors. It doesn't really matter at t=0 anyway
    // _previous_energy = _neighboring_energies[0];
};

double SpinSystem::energy() const
{
    return emap_ptr->get_config_energy(current_state);
}

std::string SpinSystem::binary_state() const
{
    int* binary_array = 0;
    binary_array = new int [params.N_spins];
    state::int_array_from_arbitrary_precision_integer_(binary_array, params.N_spins, current_state);

    std::string s = "";
    for (int ii=0; ii<params.N_spins; ii++)
    {
        s += std::to_string(binary_array[ii]);
    }

    delete[] binary_array;

    return s;
}

/* Finds the lowest energy configuration via greedy search from a config
 * through its neighbors. Returns the config index of this inherent
 structure. */
// long long SpinSystem::_help_get_inherent_structure() const
// {

//     int min_el;
//     double tmp_energy;

//     // Make a copy of the config
//     int config_copy[rtp.N_spins];
//     memcpy(config_copy, _spin_config, rtp.N_spins*sizeof(int));

//     // Neighboring energies local copy
//     double ne[rtp.N_spins];

//     // While we are not at the inherent structure, keep going
//     while (true)
//     {
//         // Compute the neighboring energies on the copy
//         _helper_fill_neighboring_energies(config_copy, rtp.N_spins, ne);

//         // If we have not reached the lowest local energy which can be reached
//         // by flipping a spin
//         min_el = min_element(ne, rtp.N_spins);
//         tmp_energy = emap_ptr->get_config_energy(
//             binary_vector_to_int(config_copy, rtp.N_spins));

//         if (ne[min_el] < tmp_energy)
//         {
//             _helper_flip_spin_(config_copy, min_el);
//         }
//         else
//         {
//             return binary_vector_to_int(config_copy, rtp.N_spins);
//         }
//     }

//     throw std::runtime_error("Unknown error in _help_get_inherent_structure");
// }


/* Updates the inherent_structure_mapping if needed and returns the config
of the inherent structure */
// long long SpinSystem::_get_inherent_structure() const
// {
//     if (rtp.memory == 0)
//     {
//         throw std::runtime_error(
//             "Inherent structure not defined for memoryless spin systems");
//     }

//     // Get the current configuration integer representations and energies
//     const long long config_int
//         = binary_vector_to_int(_spin_config, rtp.N_spins);

//     // Update the inherent structure dictionary
//     long long config_IS_int;
//     if (_ism[config_int] != -1)
//     {
//         // We've computed the inherent structure before, no need to do
//         // it again
//         config_IS_int = _ism[config_int];
//     }
//     else
//     {
//         config_IS_int = _help_get_inherent_structure();
//         _ism[config_int] = config_IS_int;
//     }
//     return config_IS_int;
// }


// void SpinSystem::_init_prev()
// {
//     prev.int_rep = _get_current_int_rep();
//     prev.energy = _get_current_energy();
//     if (rtp.memory != 0)
//     {
//         prev.int_rep_IS = _get_inherent_structure();
//         prev.energy_IS = emap_ptr->get_config_energy(prev.int_rep_IS);
//     }
// }

// void SpinSystem::_init_curr()
// {
//     curr.int_rep = _get_current_int_rep();
//     curr.energy = _get_current_energy();
//     if (rtp.memory != 0)
//     {
//         curr.int_rep_IS = _get_inherent_structure();
//         curr.energy_IS = emap_ptr->get_config_energy(curr.int_rep_IS);
//     }
// }


// double SpinSystem::get_average_neighboring_energy() const
// {
//     // Make a copy of the config if the memory != 0
//     int config_copy[rtp.N_spins];
//     double ne[rtp.N_spins];
//     if (rtp.memory != 0)
//     {
//         memcpy(config_copy, _spin_config, rtp.N_spins*sizeof(int));
//     }
//     _helper_fill_neighboring_energies(config_copy, rtp.N_spins, ne);

//     if ((rtp.memoryless_retain_last_energy == 1) and (rtp.memory == 0))
//     {
//         ne[0] = _previous_energy;
//     }

//     double tmp_e = 0.0;
//     for (int ii=0; ii<rtp.N_spins; ii++){tmp_e += ne[ii];}
//     return tmp_e / ((double) rtp.N_spins);
// }


SpinSystem::~SpinSystem()
{
    delete[] neighbors;
    delete[] neighboring_energies;
}


// Gillespie Spin system ------------------------------------------------------


// GillespieSpinSystem::GillespieSpinSystem(const RuntimeParameters rtp,
//     EnergyMapping& emap) : SpinSystem(rtp, emap)
// {
//     // Initialize the other pointers to gillespie-only required arrays
//     _delta_E = new double[rtp.N_spins];
//     _exit_rates = new double[rtp.N_spins];

//     // Initialize the normalized exit rate object
//     for (int ii=0; ii<rtp.N_spins; ii++){_normalized_exit_rates.push_back(0.0);}
// }


// double GillespieSpinSystem::_calculate_exit_rates() const
// {
//     const double current_energy = _get_current_energy();

//     for (int ii=0; ii<rtp.N_spins; ii++)
//     {
//         _delta_E[ii] = _neighboring_energies[ii] - current_energy;
//         _exit_rates[ii] = exp(-rtp.beta * _delta_E[ii]);
//         if (_exit_rates[ii] > 1.0){_exit_rates[ii] = 1.0;}
//     }
//     for (int ii=0; ii<rtp.N_spins; ii++)
//     {
//         _exit_rates[ii] = _exit_rates[ii] / ((double) rtp.N_spins);
//     }

//     double total_exit_rate = 0.0;
//     for (int ii=0; ii<rtp.N_spins; ii++)
//     {
//         total_exit_rate += _exit_rates[ii];

//     }
//     return total_exit_rate;
// }

// long double GillespieSpinSystem::step()
// {

//     const double current_energy = _get_current_energy();

//     // Update the previous state with the current information before flipping
//     _init_prev();

//     // Make a copy of the config if the memory != 0
//     int config_copy[rtp.N_spins];
//     if (rtp.memory != 0)
//     {
//         memcpy(config_copy, _spin_config, rtp.N_spins*sizeof(int));
//     }

//     // This just gets a list of random numbers if memoryless
//     _helper_fill_neighboring_energies(config_copy,
//         rtp.N_spins, _neighboring_energies);

//     // IF whatever flag is true, set one of the new neighboring energies
//     // equal to the old one
//     if ((rtp.memoryless_retain_last_energy == 1) and (rtp.memory == 0))
//     {
//         _neighboring_energies[0] = _previous_energy;
//     }

//     const double total_exit_rate = _calculate_exit_rates();

//     for (int ii=0; ii<rtp.N_spins; ii++)
//     {
//         _normalized_exit_rates[ii] = _exit_rates[ii] / total_exit_rate;
//     }

//     // Now, we make a choice of the spin to flip
//     std::discrete_distribution<int> _dist(
//         _normalized_exit_rates.begin(), _normalized_exit_rates.end());
//     const int spin_to_flip = _dist(generator);

//     // Flip that spin for real in the simulations with memory
//     if (rtp.memory != 0)
//     {
//         _flip_spin(spin_to_flip);
//     }

//     // This is the only place that the system energy and system config
//     // will be updated in the memoryless simulation!!!
//     else
//     {
//         _memoryless_system_energy = _neighboring_energies[spin_to_flip];
//         _memoryless_system_config += 1;
//     }

//     _init_curr();  // Initialize the current state

//     total_exit_rate_dist.param(
//         std::exponential_distribution<long double>::param_type(total_exit_rate));

//     _previous_energy = current_energy;

//     // Return the waiting time
//     return total_exit_rate_dist(generator) * _waiting_time_multiplier;
// }

// GillespieSpinSystem::~GillespieSpinSystem()
// {
//     delete[] _delta_E;
//     delete[] _exit_rates;
// }


// // Standard Spin system -------------------------------------------------------

// StandardSpinSystem::StandardSpinSystem(const RuntimeParameters rtp,
//     EnergyMapping& emap) : SpinSystem(rtp, emap)
// {
//     uniform_0_1_distribution.param(
//         std::uniform_real_distribution<>::param_type(0.0, 1.0));
//     spin_distribution.param(
//         std::uniform_int_distribution<>::param_type(0, rtp.N_spins - 1));
// };


// long double StandardSpinSystem::step()
// {

//     const double current_energy = _get_current_energy();

//     _init_prev();  // Initialize the current state

//     const double intermediate_energy = prev.energy;
//     const int spin_to_flip = spin_distribution(generator);

//     if (rtp.memory != 0)
//     {
//         _flip_spin(spin_to_flip);

//         // Step 4, get the proposed energy (energy of the new configuration)
//         // long long proposed_config_int = get_int_rep();
//         const double proposed_energy = _get_current_energy();

//         // Step 5, compute the difference between the energies, and find the
//         // metropolis criterion
//         const double dE = proposed_energy - intermediate_energy;
//         const double metropolis_prob = exp(-rtp.beta * dE);

//         // Step 6, sample a random number between 0 and 1.
//         const double sampled = uniform_0_1_distribution(generator);

//         // Step 7, determine whether or not to remain in this configuration or
//         // to flip back. If the randomly sampled value is less than the
//         // metropolis probability, we accept that new configuration. An easy
//         // sanity check for this is when dE is negative, then the argument of
//         // the exponent is positive and the e^(...) > 1 always; thus we always
//         // accept. However, if dE is positive, we only accept with some
//         // probability that decays exponentially quickly with the difference
//         // in energy.
//         if (sampled > metropolis_prob)  // reject
//         {
//             _flip_spin(spin_to_flip);  // flip the spin back
//         }
//     }

//     else
//     {
//         // There's no memory, so we just sample a random energy
//         const double proposed_energy = _neighboring_energies[spin_to_flip];
//         const double dE = proposed_energy - intermediate_energy;
//         const double metropolis_prob = exp(-rtp.beta * dE);
//         const double sampled = uniform_0_1_distribution(generator);
//         if (sampled <= metropolis_prob)  // accepted
//         {
//             _memoryless_system_energy = proposed_energy;
//             _memoryless_system_config += 1;

//             // Accepting a move leads to a new state and new neighbors
//             _helper_fill_neighboring_energies(_spin_config, rtp.N_spins,
//                 _neighboring_energies);

//             if (rtp.memoryless_retain_last_energy == 1)
//             {
//                 _neighboring_energies[0] = _previous_energy;
//             }
//         }
//     }

//     _previous_energy = current_energy;

//     // This is called multiple times in the loop dynamics, eventually
//     // ending with the actual value in the loop dynamics.
//     _init_curr();

//     return _waiting_time_multiplier;
// }


