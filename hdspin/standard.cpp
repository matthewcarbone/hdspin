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

#include "utils/init_utils.h"
#include "utils/general_utils.h"


// const int print_partitions = 50;



void standard(const std::string file_dump_loc, const long int N_timesteps,
    const int N_spins, const double beta, const double beta_critical,
    const int landscape)
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
    std::uniform_int_distribution<> spin_distribution(0, N_spins - 1);

    // Initialize a distribution that can pick a random number in [0, 1)
    std::uniform_real_distribution<> uniform_0_1_distribution(0.0, 1.0);

    // ========================================================================
    // Spin system & simulation ===============================================
    // ========================================================================

    // The config is an N_spins-length integer binary array indexing whether a
    // spin is up or down.
    int config[N_spins];
    initialize_spin_system(config, N_spins);

    // Initialize the energy dictionary or array for faster lookups. Note that
    // the energy array is huge and need to be explicitly allocated on the
    // heap else we will get a stackoverflow error for N ~ 20 or so.
    const int n_configs = int(pow(2, N_spins));
    double *energy_arr = new double[n_configs];
    if (landscape == 0)
    {
        initialize_energy_mapping_exponential_arr(energy_arr, N_spins,
            beta_critical);
    }
    else
    {
        initialize_energy_mapping_gaussian_arr(energy_arr, N_spins,
            beta_critical);
    }

    // The current energy is the energy of the current configuraiton BEFORE
    // stepping to the next one at the end of each algorithm step.
    double current_energy = energy_arr[binary_vector_to_int(config, N_spins)];
    double proposed_energy;

    // Terms for the inherent structure
    int config_IS_int, config_int;
    double energy_IS;

    // Initialize an array for tracking the inherent structures. This is
    // basically a mapping between the index of the array (configuration) and
    // the inherent structure configuration, the value.
    int *inherent_structure_mapping = new int[n_configs];

    // Store every entry as -1 (to indicate that none exists yet)
    for (int ii=0; ii<n_configs; ii++){inherent_structure_mapping[ii] = -1;}

    config_int = binary_vector_to_int(config, N_spins);
    config_IS_int = compute_inherent_structure(config, energy_arr, N_spins);
    inherent_structure_mapping[config_int] = config_IS_int;
    energy_IS = energy_arr[config_IS_int];


    // ========================================================================
    // Run the simulation =====================================================
    // ========================================================================

    std::ofstream outFile(file_dump_loc);

    // auto start = std::chrono::high_resolution_clock::now();
    // const double print_every = ((double) N_timesteps) / print_partitions;
    // double current_print_time = print_every;

    int spin_to_flip, sampled;
    double dE, metropolis_prob, new_energy;
    bool save_next_state = true;
    for (int timestep=0; timestep<N_timesteps; timestep++)
    {

        // Step 0: save the current results. This is more difficult for the
        // standard simulation. We only want to save when the configuration
        // changes.
        if (save_next_state == true)
        {
            outFile << timestep << " " << config_int << " " << current_energy
                << " " << config_IS_int << " " << energy_IS << "\n";
            save_next_state = false;
        }

        /*
        if (timestep > current_print_time)
        {
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration_seconds = 
                std::chrono::duration_cast<std::chrono::seconds>(stop - start);
            double duration_double_seconds = 
                std::chrono::duration<double>(duration_seconds).count();
            std::cout << "t=" << timestep << " " << "elapsed="
                << duration_double_seconds << std::endl;
            current_print_time += print_every;
        }
        */
        
        config_int = binary_vector_to_int(config, N_spins);

        // Step 2, select a random spin to flip
        spin_to_flip = spin_distribution(generator);

        // Step 3, flip that spin
        flip_spin_(config, spin_to_flip);

        // Step 4, get the proposed energy (energy of the new configuration)
        proposed_energy = energy_arr[binary_vector_to_int(config, N_spins)];

        // Step 5, compute the difference between the energies, and find the
        // metropolis criterion
        dE = proposed_energy - current_energy;
        metropolis_prob = exp(-beta * dE);

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
            new_energy = current_energy;
        }
        else  // accept, don't flip the spin back
        {
            // set the energy to the new value
            new_energy = proposed_energy;

            // Then calculate the inherent structure
            if (inherent_structure_mapping[config_int] != -1)
            {
                // We've computed the inherent structure before, no need to do
                // it again
                config_IS_int = inherent_structure_mapping[config_int];
            }
            else
            {
                config_IS_int = compute_inherent_structure(config, energy_arr,
                    N_spins);
                inherent_structure_mapping[config_int] = config_IS_int;
            }
            energy_IS = energy_arr[config_IS_int];

            // Note that if we reject then we do not save the "next" config
            save_next_state = true;
        }

        current_energy = new_energy;
    }

    delete[] energy_arr;
    delete[] inherent_structure_mapping;
}
