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

    // Default assume rejected
    bool accepted = false;

    // Distributions are cheap, so we can safely re-initialize them at
    // every step. The generator, however, has a lot of overhead and should
    // be passed.
    // Randomly pick a random number in [0, 1)
    std::uniform_real_distribution<> uniform_0_1_distribution(0.0, 1.0);

    int spin_to_flip = 0;
    while (spin_to_flip < rtp.N_spins)
    {
        init_prev_();  // Initialize the current state

        // Step 3, flip that spin
        if (rtp.loop_dynamics == 1)
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
        long double proposed_energy = get_current_energy();

        // Step 5, compute the difference between the energies, and find the
        // metropolis criterion
        long double dE = proposed_energy - prev.energy;
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

            /*
            // Step the standard config psi counter
            psi_config_counter.step(sys.waiting_time, false);
            sys.waiting_time = 0.0;

            // Update the inherent structure values
            inh.x_prev = inh.x;
            inh.e_prev = inh.e;
            inh.x = query_inherent_structure(params.N_spins, config,
                energy_arr, inherent_structure_mapping);
            inh.e = energy_arr[inh.x];

            // This is a tricky update for the inherent structure, since it
            // will have a different waiting time than the normal
            // configuration, as it may not change even though the normal
            // configuration does.
            if (inh.x == inh.x_prev){inh.waiting_time += 1.0;}
            else
            {
                // Step the inherent structure psi config counter
                psi_config_counter.step(inh.waiting_time, true);
                inh.waiting_time = 0.0;
            }
            */

            accepted = true;
        }

        init_curr_();

        if (rtp.loop_dynamics == 1){break;}

        spin_to_flip++;
    }

    // If no change is ever accepted, increment
    if (!accepted)
    {
        time_in_config += 1.0;
        time_in_config_IS += 1.0;
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
