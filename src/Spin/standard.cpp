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
    SpinSystem(rtp) {};


bool StandardSpinSystem::step_()
{
    // The loop_dynamics parameter is really important here. If it is equal
    // to 1, we have loop dynamics, where every spin is checked in order. If
    // it is not equal to 1, we have either standard (0) or div-N (2) dynamics,
    // in which case only one spin is checked randomly. The difference in
    // waiting time is handled outside of this module, in sim.cpp.

    // Default assume rejected
    bool accepted = false;

    // Distributions are cheap, so we can safely re-initialize them at
    // every step. The generator, however, has a lot of overhead and should
    // be passed.
    // Randomly pick a random number in [0, 1)
    std::uniform_real_distribution<> uniform_0_1_distribution(0.0, 1.0);

    // We assume in the outer loop that we're running loop dynamics ...
    int spin_to_flip = 0;

    // We need to differentiate between the loop dynamics and the standard ones
    // when initializing the previous spin. In the loop dynamics, the previous
    // energy should only be recorded once at the beginning.
    init_prev_();  // Initialize the current state

    double intermediate_energy = prev.energy;

    while (spin_to_flip < rtp.N_spins)
    {
        // ... but if we're not, randomly sample the spin, overriding the
        // for loop
        if (rtp.loop_dynamics != 1)
        {
            // Randomly pick a spin from 0 -> N - 1
            std::uniform_int_distribution<> spin_distribution(
                0, rtp.N_spins - 1);
            // Select a random spin to flip (override the outer loop)
            spin_to_flip = spin_distribution(generator);
        }
        flip_spin_(spin_to_flip);

        // Step 4, get the proposed energy (energy of the new configuration)
        // long long proposed_config_int = get_int_rep();
        double proposed_energy = get_current_energy();

        // Step 5, compute the difference between the energies, and find the
        // metropolis criterion
        double dE = proposed_energy - intermediate_energy;
        double metropolis_prob = exp(-rtp.beta * dE);

        // Step 6, sample a random number between 0 and 1.
        double sampled = uniform_0_1_distribution(generator);

        // Step 7, determine whether or not to remain in this configuration or
        // to flip back. If the randomly sampled value is less than the
        // metropolis probability, we accept that new configuration. An easy
        // sanity check for this is when dE is negative, then the argument of
        // the exponent is positive and the e^(...) > 1 always; thus we always
        // accept. However, if dE is positive, we only accept with some
        // probability that decays exponentially quickly with the difference
        // in energy.
        if (sampled > metropolis_prob)  // reject
        {
            flip_spin_(spin_to_flip);  // flip the spin back
        }

        else  // accept, don't flip the spin back
        {
            accepted = true;
        }

        // This is called multiple times in the loop dynamics, eventually
        // ending with the actual value in the loop dynamics.
        init_curr_();
        intermediate_energy = curr.energy;

        // Once again, if we're not running the loop over N dynamics, break
        // here.
        if (rtp.loop_dynamics != 1){break;}

        spin_to_flip++;
    }

    // If no change is ever accepted, increment the time in config counters
    if (!accepted)
    {
        // In the case of standard dynamics, and standard loop dynamics, this
        // counter increments by 1
        if (rtp.loop_dynamics != 2)
        {
            time_in_config += 1.0;
            time_in_config_IS += 1.0;
        }

        // In the case of the div-N dynamics, this means rejection has only
        // occurred for a single div-N timestep, which is of course 1/N.
        else
        {
            time_in_config += 1.0 / rtp.N_spins;
            time_in_config_IS += 1.0 / rtp.N_spins;
        }
    }

    // Else, increment n_accept
    else
    {
        n_accept += 1;

        // Check for possible (although unlikely) overflow
        assert(n_accept > 0);
    }

    return accepted;
}
