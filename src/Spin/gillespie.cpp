/* Core system utilities.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <random>

#include "Utils/structures.h"
#include "Spin/base.h"
#include "Spin/gillespie.h"
#include "Utils/utils.h"


GillespieSpinSystem::GillespieSpinSystem(const RuntimeParameters rtp) :
    SpinSystem(rtp)
{
    // Initialize the other pointers to gillespie-only required arrays
    neighboring_energies = new double[rtp.N_spins];
    delta_E = new double[rtp.N_spins];
    exit_rates = new double[rtp.N_spins];

    // Initialize the normalized exit rate object
    for (int ii=0; ii<rtp.N_spins; ii++){normalized_exit_rates.push_back(0.0);}
}


double GillespieSpinSystem::_calculate_exit_rates() const
{
    double current_energy = _get_current_energy();
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        delta_E[ii] = neighboring_energies[ii] - current_energy;
        exit_rates[ii] = exp(-rtp.beta * delta_E[ii]);
        if (exit_rates[ii] > 1.0){exit_rates[ii] = 1.0;}
    }
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        exit_rates[ii] = exit_rates[ii] / ((double) rtp.N_spins);
    }

    double total_exit_rate = 0.0;
    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        total_exit_rate += exit_rates[ii];

    }
    return total_exit_rate;
}

long double GillespieSpinSystem::step_()
{
    // Update the previous state with the current information before flipping
    _init_prev();

    _calculate_neighboring_energies();
    const double total_exit_rate = _calculate_exit_rates();

    for (int ii=0; ii<rtp.N_spins; ii++)
    {
        normalized_exit_rates[ii] = exit_rates[ii] / total_exit_rate;
    }

    // Now, we make a choice of the spin to flip
    std::discrete_distribution<int> _dist(
        normalized_exit_rates.begin(), normalized_exit_rates.end());
    const int spin_to_flip = _dist(generator);

    _flip_spin_(spin_to_flip);
    if (rtp.memoryless == 1){spin_config_energy = _get_random_energy();}

    _init_curr();  // Initialize the current state

    n_accept += 1;  // Every step is accepted using Gillespie

    total_exit_rate_dist.param(std::exponential_distribution<long double>::param_type(total_exit_rate));

    return total_exit_rate_dist(generator) * _waiting_time_multiplier;
}

GillespieSpinSystem::~GillespieSpinSystem()
{
    delete[] delta_E;
    delete[] exit_rates;
    delete[] neighboring_energies;
}
