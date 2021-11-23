/* Core system utilities.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <cassert>
#include <random>

#include "Spin/standard.h"
#include "Utils/structures.h"

StandardSpinSystem::StandardSpinSystem(const RuntimeParameters rtp) :
    SpinSystem(rtp)
{
    uniform_0_1_distribution.param(std::uniform_real_distribution<>::param_type(0.0, 1.0));
    spin_distribution.param(std::uniform_int_distribution<>::param_type(0, rtp.N_spins - 1));
};


long double StandardSpinSystem::step_()
{
    // The loop_dynamics parameter is really important here. If it is equal
    // to 1, we have loop dynamics, where every spin is checked in order. If
    // it is not equal to 1, we have either standard (0) or div-N (2) dynamics,
    // in which case only one spin is checked randomly. The difference in
    // waiting time is handled outside of this module, in sim.cpp.

    // We need to differentiate between the loop dynamics and the standard ones
    // when initializing the previous spin. In the loop dynamics, the previous
    // energy should only be recorded once at the beginning.

    _init_prev();  // Initialize the current state

    const double intermediate_energy = prev.energy;

    const int spin_to_flip = spin_distribution(generator);
    _flip_spin_(spin_to_flip);

    // Step 4, get the proposed energy (energy of the new configuration)
    // long long proposed_config_int = get_int_rep();
    const double proposed_energy = _get_current_energy();

    // Step 5, compute the difference between the energies, and find the
    // metropolis criterion
    const double dE = proposed_energy - intermediate_energy;
    const double metropolis_prob = exp(-rtp.beta * dE);

    // Step 6, sample a random number between 0 and 1.
    const double sampled = uniform_0_1_distribution(generator);

    // Step 7, determine whether or not to remain in this configuration or
    // to flip back. If the randomly sampled value is less than the
    // metropolis probability, we accept that new configuration. An easy
    // sanity check for this is when dE is negative, then the argument of
    // the exponent is positive and the e^(...) > 1 always; thus we always
    // accept. However, if dE is positive, we only accept with some
    // probability that decays exponentially quickly with the difference
    // in energy.
    bool accepted;
    if (sampled > metropolis_prob)  // reject
    {
        _flip_spin_(spin_to_flip);  // flip the spin back
        accepted = false;
    }

    else  // accept, don't flip the spin back
    {
        accepted = true;
    }

    // This is called multiple times in the loop dynamics, eventually
    // ending with the actual value in the loop dynamics.
    _init_curr();

    // If no change is ever accepted, increment the time in config counters
    if (!accepted)
    {
        time_in_config += _waiting_time_multiplier;
        time_in_config_IS += _waiting_time_multiplier;
    }

    // Else, increment n_accept
    else
    {
        n_accept += 1;

        // Check for possible (although unlikely) overflow
        assert(n_accept > 0);
    }

    return _waiting_time_multiplier;
}
