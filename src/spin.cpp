#include <unistd.h>
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

    unsigned int* spin_config = 0;
    spin_config = new unsigned int [params.N_spins];

    std::bernoulli_distribution _bernoulli_distribution;

    for (unsigned int ii=0; ii<params.N_spins; ii++)
    {
        spin_config[ii] = _bernoulli_distribution(generator);
    }

    utils::arbitrary_precision_integer_from_int_array_(spin_config, params.N_spins, current_state);

    delete[] spin_config;

    if (params.dynamics == "standard"){_init_standard();}
    else if (params.dynamics == "gillespie"){_init_gillespie();}
    else
    {
        throw std::runtime_error("Uknown dynamics during setup");
    }

}

void SpinSystem::_init_previous_state_()
{
    _prev.state = current_state;
    _prev.energy = energy();
}

std::string SpinSystem::get_previous_state_string_rep() const
{
    return utils::string_rep_from_arbitrary_precision_integer(_prev.state, params.N_spins);
}


void SpinSystem::_init_current_state_()
{
    _curr.state = current_state;
    _curr.energy = energy();
}

std::string SpinSystem::get_current_state_string_rep() const
{
    return utils::string_rep_from_arbitrary_precision_integer(_curr.state, params.N_spins);
}

SpinSystem::SpinSystem(const utils::SimulationParameters params,
    EnergyMapping& emap) : params(params)
{
    emap_ptr = &emap;
    _first_time_state_initialization_();
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
    return utils::string_rep_from_arbitrary_precision_integer(current_state, params.N_spins);
}


void SpinSystem::_init_gillespie()
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

void SpinSystem::_teardown_gillespie()
{
    delete[] _exit_rates;
    delete[] _neighbors;
    delete[] _neighboring_energies;
}


double SpinSystem::_calculate_exit_rates(const double current_energy) const
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

double SpinSystem::_step_gillespie()
{

    // Initialize the current state as _prev
    _init_previous_state_();

    // Get the current energy of the state
    const double current_energy = _prev.energy;

    // Get the neighboring states
    utils::get_neighbors_(_neighbors, current_state, params.N_spins);

    // Populate the neighboring energies
    emap_ptr->get_config_energies_array_(_neighbors, _neighboring_energies, params.N_spins);

    const double total_exit_rate = _calculate_exit_rates(current_energy);

    for (unsigned int ii=0; ii<params.N_spins; ii++)
    {
        _normalized_exit_rates[ii] = _exit_rates[ii] / total_exit_rate;
    }

    // Now, we make a choice of the spin to flip
    std::discrete_distribution<int> _dist(
        _normalized_exit_rates.begin(), _normalized_exit_rates.end());

    // The spin to flip is actually on the "opposite side" because of how
    // bits work
    const unsigned int spin_to_flip = _dist(generator);

    // And always flip that spin in a Gillespie simulation
    current_state = utils::flip_bit(current_state, spin_to_flip, params.N_spins);

    // Initialize the current state
    _init_current_state_();

    // Calculate the waiting time
    total_exit_rate_dist.param(
        std::exponential_distribution<double>::param_type(total_exit_rate));

    // Return the waiting time which is generally != 1
    sim_stats.acceptances += 1;  // Gillespie always accepts! =)
    return total_exit_rate_dist(generator);
}


void SpinSystem::_init_standard()
{
    uniform_0_1_distribution.param(
        std::uniform_real_distribution<>::param_type(0.0, 1.0));
    spin_distribution.param(
        std::uniform_int_distribution<>::param_type(0, params.N_spins - 1));
}

double SpinSystem::_step_standard()
{

    // Initialize the current state as _prev
    _init_previous_state_();

    // Get the current energy of the state
    const double current_energy = _prev.energy;

    // Select a random spin to flip
    const unsigned int bit_to_flip = spin_distribution(generator);

    // Flip the current state into its new one
    const ap_uint<PRECISON> possible_state = utils::flip_bit(current_state, bit_to_flip, params.N_spins);

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
        sim_stats.acceptances += 1;
    }
    else
    {
        sim_stats.rejections += 1;
    }

    _init_current_state_();

    return 1.0;
}

void SpinSystem::_teardown_standard(){;}


void SpinSystem::summarize()
{
    printf("Acceptances/rejections: %lli/%lli\n", sim_stats.acceptances, sim_stats.rejections);
}

double SpinSystem::step()
{
    auto t_start = std::chrono::high_resolution_clock::now();    
    double waiting_time;
    if (params.dynamics == "standard")
    {
        waiting_time = _step_standard();
    }
    else if (params.dynamics == "gillespie")
    {
        waiting_time = _step_gillespie();
    }
    else
    {
        throw std::runtime_error("Uknown dynamics during step");
    }
    const double duration = utils::get_time_delta(t_start);
    sim_stats.total_wall_time += duration;
    sim_stats.total_waiting_time += waiting_time;
    sim_stats.total_steps += 1;
    return waiting_time;
}

SpinSystem::~SpinSystem()
{
    if (params.dynamics == "standard"){_teardown_standard();}
    else if (params.dynamics == "gillespie"){_teardown_gillespie();}

    // Else both to be safe
    else{_teardown_standard(); _teardown_gillespie();}
}
