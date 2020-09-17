/* Core local spin system algorithm.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <math.h>
#include <string>
#include <fstream>      // std::ofstream
#include <assert.h>
#include <stdio.h>

#include "utils/init_utils.h"
#include "utils/general_utils.h"
#include "utils/grid_utils.h"
#include "utils/structure_utils.h"


void gillespie(const FileNames fnames, const RuntimeParameters params)
{

    // ========================================================================
    // Random =================================================================
    // ========================================================================

    // Initialize the MT random number generator and seed with random_device
    std::mt19937 generator;
    unsigned int seed = std::random_device{}();
    generator.seed(seed);

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
    const long long N_timesteps = ipow(10, params.log_N_timesteps);
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

    // Vector of the neighboring energies which is rewritten at every step of
    // the while loop. Also a vector of the dE values, exit rates...
    double neighboring_energies[params.N_spins];
    double exit_rates[params.N_spins];
    double delta_E[params.N_spins];
    double total_exit_rate;

    // Initialize an array for tracking the inherent structures. This is
    // basically a mapping between the index of the array (configuration) and
    // the inherent structure configuration, the value.
    long long *inherent_structure_mapping = new long long[n_configs];

    // Store every entry as -1 (to indicate that none exists yet)
    for (long long ii=0; ii<n_configs; ii++){
        inherent_structure_mapping[ii] = -1;
    }

    // ========================================================================
    // Run the simulation =====================================================
    // ========================================================================

    long double current_time = 0.0;
    long double waiting_time;

    // The config index will iterate very time we change configurations, and
    // thus it does not represent the configuration itself, but pretends that
    // every configuration, even if it is revisited, is different.
    long long n_accept = 0;

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

    SystemInformation sys, inh;

    // Initialize the original values
    sys.x = binary_vector_to_int(config, params.N_spins);
    sys.x_prev = sys.x;
    sys.e = energy_arr[sys.x];
    sys.e_prev = sys.e;

    // Do the same for the inherent structure
    inh.x = query_inherent_structure(params.N_spins, config, energy_arr,
        inherent_structure_mapping);
    inh.x_prev = inh.x;
    inh.e = energy_arr[inh.x];
    inh.e_prev = inh.e;

    while (true)
    {

        // --------------------------------------------------------------------
        // ---------------------------- ENGINE --------------------------------
        // --------------------------------------------------------------------

        // Step 1: get the neighboring energies by filling the relevant object
        get_neighboring_energies(config, energy_arr, neighboring_energies,
            params.N_spins);

        // Step 2: get the exit rates and dE values
        get_exit_rates(sys.e, params.beta, neighboring_energies,
            exit_rates,  delta_E, params.N_spins);
        total_exit_rate = 0.0;
        for (int ii=0; ii<params.N_spins; ii++)
        {
            total_exit_rate += exit_rates[ii];
        }

        // Step 3: initialize an exponential distribution to sample the
        // waiting time from
        std::exponential_distribution<long double>
            tmp_exp_dist(total_exit_rate);
        waiting_time = tmp_exp_dist(generator);

        // At this point: we understand that the tracer is in the configuration
        // it was in at the beginning of this step, and that it is in that
        // configuration for `waiting_time` time, and the new current_time
        // after stepping is `current_time` + `waiting_time`.

        // Step 4: update the current time of the simulation clock
        current_time += waiting_time;

        // Step 6: step to the next state and store the proposed (new) energy
        step_next_state_(config, exit_rates, total_exit_rate, params.N_spins,
            generator);

        //         -------------------------------------------------
        //         ---------------- STEP TRACKERS ------------------
        //         -------------------------------------------------

        // Append trackers. Note that Gillespie dynamics are different than
        // standard in the order in which we update the grids, so the grids are
        // actually stepped before the sys/inh objects are updated.  
        energy_grid.step(current_time, sys, inh);
        aging_config_grid.step(current_time, n_accept, sys.x, inh.x);
        psi_config_counter.step(waiting_time, false);  // Step standard

        // This is a tricky update for the inherent structure, since it
        // will have a different waiting time than the normal
        // configuration, as it may not change even though the normal
        // configuration does.
        if (inh.x == inh.x_prev){inh.waiting_time += waiting_time;}
        else
        {
            // Step the inherent structure psi config counter
            psi_config_counter.step(inh.waiting_time, true);  // Step IS
            inh.waiting_time = 0.0;
        }

        //         -------------------------------------------------
        //         -------------- DONE STEP TRACKERS ---------------
        //         -------------------------------------------------

        // Update values for the system
        sys.x_prev = sys.x;
        sys.x = binary_vector_to_int(config, params.N_spins);
        sys.e_prev = sys.e;
        sys.e = energy_arr[sys.x];
        

        // Update the inherent structure values
        inh.x_prev = inh.x;
        inh.e_prev = inh.e;
        inh.x = query_inherent_structure(params.N_spins, config,
            energy_arr, inherent_structure_mapping);
        inh.e = energy_arr[inh.x];

        

        // --------------------------------------------------------------------
        // ----------------------- ENGINE FINISH ------------------------------
        // --------------------------------------------------------------------

        n_accept += 1;

        // Check for possible (although unlikely) overflow
        assert(current_time > 0.0);
        assert(n_accept > 0);

        if (current_time >= N_timesteps){break;}

    }

    // Close the outfiles and write to disk when not doing so dynamically
    energy_grid.close_outfile();
    psi_config_counter.write_to_disk(fnames.psi_config);
    aging_config_grid.close_outfile();

    delete[] energy_arr;
    delete[] inherent_structure_mapping;
}
