/* Core local spin system algorithm.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <math.h>
#include <string>
#include <fstream>      // std::ofstream
#include <assert.h>
#include <chrono>
#include <stdio.h>

#include "utils/init_utils.h"
#include "utils/general_utils.h"
#include "utils/grid_utils.h"
#include "utils/structure_utils.h"


struct SystemInformation
{
    long long x, x_prev;
    long double e, e_prev;

    // The amount of time the system has spent in some configuration
    long double waiting_time = 0.0;
};


void standard(const FileNames fnames, const RuntimeParameters params)
{

    /* We first initialize counters, trackers and grids for the various
    observables we are interested in recording.*/

    // ========================================================================
    // Random =================================================================
    // ========================================================================

    // Initialize the MT random number generator and seed with random_device
    std::mt19937 generator;
    unsigned int seed = std::random_device{}();
    generator.seed(seed);

    // Initialize a distribution that can randomly pick a spin from 0 -> N - 1
    std::uniform_int_distribution<> spin_distribution(0, params.N_spins - 1);

    // Initialize a distribution that can pick a random number in [0, 1)
    std::uniform_real_distribution<> uniform_0_1_distribution(0.0, 1.0);

    // ========================================================================
    // Spin system & simulation ===============================================
    // ========================================================================

    // The config is an N_spins-length integer binary array indexing whether a
    // spin is up or down.
    int config[params.N_spins];
    initialize_spin_system(config, params.N_spins);

    // Initialize the energy dictionary or array for faster lookups. Note that
    // the energy array is huge and need to be explicitly allocated on the
    // heap else we will get a stackoverflow error for N ~ 20 or so.
    const long long n_configs = ipow(2, params.N_spins);
    const long long N_timesteps = ipow(10, params.log_N_timesteps) + 2;

    double *energy_arr = new double[n_configs];
    if (params.landscape == 0)
    {
        initialize_energy_mapping_exponential_arr(energy_arr, n_configs,
            params.beta_critical);
    }
    else
    {
        initialize_energy_mapping_gaussian_arr(energy_arr, n_configs,
            params.N_spins, params.beta_critical);
    }

    // The current energy is the energy of the current configuraiton BEFORE
    // stepping to the next one at the end of each algorithm step.
    double proposed_energy;
    long long proposed_config_int;

    // Initialize an array for tracking the inherent structures. This is
    // basically a mapping between the index of the array (configuration) and
    // the inherent structure configuration, the value.
    long long *inherent_structure_mapping = new long long[n_configs];

    // Store every entry as -1 (to indicate that none exists yet)
    for (long long ii=0; ii<n_configs; ii++)
    {
        inherent_structure_mapping[ii] = -1;
    }


    // ========================================================================
    // Run the simulation =====================================================
    // ========================================================================

    int spin_to_flip;
    double dE, metropolis_prob, sampled;

    // Define all the observable trackers ---------------------------------
    // Energy
    EnergyGrid energy_grid(fnames.grids_directory);
    energy_grid.open_outfile(fnames.energy);
    // Psi config
    PsiConfigCounter psi_config_counter(params.log_N_timesteps);
    // Pi/Aging config
    AgingConfigGrid aging_config_grid(fnames.grids_directory);
    aging_config_grid.open_outfile(fnames.aging_config_1,
        fnames.aging_config_2);

    // Initialize the original values
    SystemInformation sys, sys_IS;
    sys.x = binary_vector_to_int(config, params.N_spins);
    sys.x_prev = sys.x;
    sys.e = energy_arr[sys.x];
    sys.e_prev = sys.e;

    // Do the same for the inherent structure
    sys_IS.x = query_inherent_structure(params.N_spins, config,
        energy_arr, inherent_structure_mapping);
    sys_IS.x_prev = sys_IS.x;
    sys_IS.e = energy_arr[sys_IS.x];
    sys_IS.e_prev = sys_IS.e;

    // This index will iterate very time we accept a new configuration, and
    // thus it does not represent the configuration itself, but pretends that
    // every configuration, even if it is revisited, is different.
    long long n_accepted = 0;

    for (long long timestep=0; timestep<N_timesteps; timestep++)
    {
        
        // --------------------------------------------------------------------
        // ---------------------------- ENGINE --------------------------------
        // --------------------------------------------------------------------

        // Step 2, select a random spin to flip
        spin_to_flip = spin_distribution(generator);

        // Step 3, flip that spin
        flip_spin_(config, spin_to_flip);

        // Step 4, get the proposed energy (energy of the new configuration)
        proposed_config_int = binary_vector_to_int(config, params.N_spins);
        proposed_energy = energy_arr[proposed_config_int];

        // Step 5, compute the difference between the energies, and find the
        // metropolis criterion
        dE = proposed_energy - sys.e;
        metropolis_prob = exp(-params.beta * dE);

        // Step 6, sample a random number between 0 and 1.
        sampled = uniform_0_1_distribution(generator);

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
            flip_spin_(config, spin_to_flip);  // flip the spin back
            sys.waiting_time += 1.0;
            sys_IS.waiting_time += 1.0;
        }
        else  // accept, don't flip the spin back
        {
            // Update the values for the system and inherent structure
            sys.x_prev = sys.x;
            sys.x = proposed_config_int;
            sys.e_prev = sys.e;
            sys.e = proposed_energy;

            psi_config_counter.step(sys.waiting_time);
            sys.waiting_time = 0.0;

            // Update the inherent structure values
            sys_IS.x_prev = sys_IS.x;
            sys_IS.e_prev = sys_IS.e;
            sys_IS.x = query_inherent_structure(params.N_spins, config,
                energy_arr, inherent_structure_mapping);
            sys_IS.e = energy_arr[sys_IS.x];

            // This is a tricky update for the inherent structure, since it
            // will have a different waiting time than the normal
            // configuration, as it may not change even though the normal
            // configuration does.
            if (sys_IS.x == sys_IS.x_prev){sys_IS.waiting_time += 1.0;}
            else
            {
                // NEED TO APPEND A NEW COUNTER HERE FOR PSI IS!!!
                sys_IS.waiting_time = 0.0;
            }

            n_accepted += 1;

            // Check for possible (although unlikely) overflow
            assert(n_accepted > 0);
        }

        // --------------------------------------------------------------------
        // ----------------------- ENGINE FINISH ------------------------------
        // --------------------------------------------------------------------

        // Append trackers
        energy_grid.step(timestep, sys.x, sys_IS.x, sys.e, sys_IS.e);
        aging_config_grid.step(timestep, n_accepted, sys.x, sys_IS.x);

    }

    // Close the outfiles and write to disk when not doing so dynamically
    energy_grid.close_outfile();
    psi_config_counter.write_to_disk(fnames.psi_config);
    aging_config_grid.close_outfile();

    delete[] energy_arr;
    delete[] inherent_structure_mapping;
}
