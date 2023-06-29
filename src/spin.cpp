#include <cstring>
#include <random>

#include "utils.h"
#include "spin.h"


void SpinSystem::_first_time_state_initialization_()
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

// void SpinSystem::_fill_neighbors_()
// {
//     state::get_neighbors_(neighbors, current_state, params.N_spins);
// }

// void SpinSystem::_fill_neighboring_energies_()
// {
//     for (int ii=0; ii<params.N_spins; ii++)
//     {
//         neighboring_energies[ii] = emap_ptr->get_config_energy(neighbors[ii]);
//     }
// }

void SpinSystem::_init_previous_state_()
{
    _prev.state = current_state;
    _prev.energy = energy();
}

void SpinSystem::_init_current_state_()
{
    _curr.state = current_state;
    _curr.energy = energy();
}

SpinSystem::SpinSystem(const parameters::SimulationParameters params,
    EnergyMapping& emap) : params(params)
{
    emap_ptr = &emap;

    _first_time_state_initialization_();

    // std::uniform_int_distribution<> int_dist(0, params.N_spins - 1);

    // if (params.divN == 1){_waiting_time_multiplier = 1.0 / params.N_spins;}


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

void SpinSystem::set_state(ap_uint<PRECISON> state)
{
    current_state = state;
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


// SpinSystem::~SpinSystem()
// {
//     delete[] neighbors;
//     delete[] neighboring_energies;
// }


// Gillespie Spin system ------------------------------------------------------

GillespieSpinSystem::GillespieSpinSystem(const parameters::SimulationParameters params,
    EnergyMapping& emap) : SpinSystem(params, emap)
{
    // Initialize the other pointers to gillespie-only required arrays
    _exit_rates = new double[params.N_spins];
    _neighbors = new ap_uint<PRECISON>[params.N_spins];
    _neighboring_energies = new double[params.N_spins];

    // Initialize the normalized exit rate object
    for (int ii=0; ii<params.N_spins; ii++)
    {
        _normalized_exit_rates.push_back(0.0);
    }
}


double GillespieSpinSystem::_calculate_exit_rates(const double current_energy) const
{
    double dE;
    for (int ii=0; ii<params.N_spins; ii++)
    {
        dE = _neighboring_energies[ii] - current_energy;
        _exit_rates[ii] = exp(-params.beta * dE);
        if (_exit_rates[ii] > 1.0){_exit_rates[ii] = 1.0;}
    }
    for (int ii=0; ii<params.N_spins; ii++)
    {
        _exit_rates[ii] = _exit_rates[ii] / ((double) params.N_spins);
    }

    double total_exit_rate = 0.0;
    for (int ii=0; ii<params.N_spins; ii++)
    {
        total_exit_rate += _exit_rates[ii];

    }
    return total_exit_rate;
}

long double GillespieSpinSystem::step()
{

    // Initialize the current state as _prev
    _init_previous_state_();

    // Get the current energy of the state
    const double current_energy = _prev.energy;

    // Get the neighboring states
    state::get_neighbors_(_neighbors, current_state, params.N_spins);

    // Populate the neighboring energies
    emap_ptr->get_config_energies_array(_neighbors, _neighboring_energies, params.N_spins);

    const double total_exit_rate = _calculate_exit_rates(current_energy);

    for (int ii=0; ii<params.N_spins; ii++)
    {
        _normalized_exit_rates[ii] = _exit_rates[ii] / total_exit_rate;
    }

    // Now, we make a choice of the spin to flip
    std::discrete_distribution<int> _dist(
        _normalized_exit_rates.begin(), _normalized_exit_rates.end());
    const int spin_to_flip = params.N_spins - 1 - _dist(generator);
    // std::cout << "flipping: " << spin_to_flip << std::endl;

    // And always flip that spin in a Gillespie simulation
    // const ap_uint<PRECISON> c1 = current_state;
    // const std::string s1 = binary_state();
    current_state = state::flip_bit(current_state, spin_to_flip, params.N_spins);
    // const ap_uint<PRECISON> c2 = current_state;
    // const std::string s2 = binary_state();
    // for (int ii=0; ii<params.N_spins; ii++){
    //     std::cout << _normalized_exit_rates[ii] << std::endl;
    // }
    // std::cout << c1 << " " << c2 << std::endl;
    // std::cout << s1 << "->" << s2 << std::endl;

    // Initialize the current state
    _init_current_state_();

    // Calculate the waiting time
    total_exit_rate_dist.param(
        std::exponential_distribution<long double>::param_type(total_exit_rate));

    // Return the waiting time which is generally != 1
    return total_exit_rate_dist(generator);
}

GillespieSpinSystem::~GillespieSpinSystem()
{
    delete[] _exit_rates;
    delete[] _neighbors;
    delete[] _neighboring_energies;
}


// Standard Spin system -------------------------------------------------------

StandardSpinSystem::StandardSpinSystem(const parameters::SimulationParameters params, EnergyMapping& emap) : SpinSystem(params, emap)
{
    uniform_0_1_distribution.param(
        std::uniform_real_distribution<>::param_type(0.0, 1.0));
    spin_distribution.param(
        std::uniform_int_distribution<>::param_type(0, params.N_spins - 1));
};


long double StandardSpinSystem::step()
{

    // Initialize the current state as _prev
    _init_previous_state_();

    // Get the current energy of the state
    const double current_energy = _prev.energy;

    // Select a random spin to flip
    const int spin_to_flip = spin_distribution(generator);

    // Flip the current state into its new one
    const ap_uint<PRECISON> possible_state = state::flip_bit(current_state, spin_to_flip, params.N_spins);

    // Get the proposed energy (energy of the new configuration)
    const double proposed_energy = emap_ptr->get_config_energy(possible_state);

    // Compute the difference between the energies, and find the metropolis
    // selection criterion
    const double dE = proposed_energy - current_energy;
    const double metropolis_prob = exp(-params.beta * dE);

    // Sample a random number between 0 and 1
    const double sampled = uniform_0_1_distribution(generator);

    // Determine whether or not to remain in this configuration or
    // to flip back. If the randomly sampled value is less than the
    // metropolis probability, we accept that new configuration. An easy
    // sanity check for this is when dE is negative, then the argument of
    // the exponent is positive and the e^(...) > 1 always; thus we always
    // accept. However, if dE is positive, we only accept with some
    // probability that decays exponentially quickly with the difference
    // in energy.
    if (sampled <= metropolis_prob)  // accept
    {
        current_state = possible_state;
        counter.acceptances += 1;
    }
    else
    {
        counter.rejections += 1;
    }

    _init_current_state_();

    return 1.0;
}


void StandardSpinSystem::summarize()
{
    printf("Acceptances/rejections: %lli/%lli\n", counter.acceptances, counter.rejections);
}

