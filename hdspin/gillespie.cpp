/* Core local spin system algorithm.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <math.h>
#include <string>

#include "utils/grid_utils.h"
#include "utils/init_utils.h"
#include "utils/general_utils.h"
#include "utils/file_utils.h"


void append_energy_tracker_(const double current_energy,
    const double current_time, const double *time_grid,
    double *energy_recorder, int *energy_pointer, const int ngrid)
{
    if (*energy_pointer >= ngrid){return;}
    if (current_time < time_grid[*energy_pointer]){return;}

    while (current_time >= time_grid[*energy_pointer])
    {
        energy_recorder[*energy_pointer] = current_energy;
        *energy_pointer += 1;

        if (*energy_pointer >= ngrid){return;}
    }
}


void append_psi_recorder_(int *counter, const double time,
    const int len_psi_grid)
{

    int idx;
    if (time < 1e-10){idx = 0;}
    else
    {
        idx = floor(log2(time));
        if(idx < 0){idx = 0;}
    }

    // In the case that we encounter a situation where the tracer exceeded the
    // limit, we simply do not record it, as these are fringe cases.
    if (idx >= len_psi_grid){return;}

    counter[idx] += 1;
}


int append_Pi_grid_recorder_(const double current_time,
    const double current_energy, const double threshold_energy,
    const int tracker_index, int tracking_pointer, const int total_gridpoints,
    const double *time_pi_grid, int *recorder_pi_grid)
{

    // Early break conditions; we don't want to accidentally crash by accessing
    // elements of the array that don't exist.
    if (tracking_pointer > total_gridpoints - 1)
    {
        return tracking_pointer;
    }

    while (current_time >= time_pi_grid[tracking_pointer])
    {
        recorder_pi_grid[tracking_pointer] = tracker_index;
        if (current_energy < threshold_energy)
        {
            // Get's a negative sign if in a basin
            recorder_pi_grid[tracking_pointer] = 
                -recorder_pi_grid[tracking_pointer];
        }

        tracking_pointer += 1;

        if (tracking_pointer > total_gridpoints - 1)
        {
            return tracking_pointer;
        }
    }

    return tracking_pointer;
}


void gillespie(const int index, const std::string file_dump_loc,
    const long int N_timesteps, const int N_spins, const int n_samp,
    const double beta, const double beta_critical, const int landscape,
    const double thresh_S, const double thresh_E, const int print_every,
    const double DW)
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
    double current_time = 0.0;

    // ========================================================================
    // Energy and related =====================================================
    // ========================================================================

    // Indexes the next-to-be-recorded position on the time_energy_grid. In
    // other words, time_energy_grid[energy_pointer] will be some float
    // representing the time at which to record the energy next on the grid.
    // This also means that time_energy_grid[energy_pointer - 1] will have
    // already been recorded.
    int energy_pointer = 0;

    // Initialize the grid which holds the time-values at which to record the
    // energies. This grid is logarithmically spaced in time.
    double time_energy_grid[n_samp];
    const double base = 10.0;
    fill_pyLogspace(time_energy_grid, 0, log10(N_timesteps), n_samp, base);

    // Initialize the recorder for the energy
    double recorder_energy_grid[n_samp];

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

    // Vector of the neighboring energies which is rewritten at every step of
    // the while loop. Also a vector of the dE values, exit rates...
    double neighboring_energies[N_spins];
    double exit_rates[N_spins];
    double delta_E[N_spins];
    double total_exit_rate;
    double waiting_time;

    // ========================================================================
    // Psi ====================================================================
    // ========================================================================

    // The psi counters log how long the tracers spend in a single basin, or
    // in a single configuration.
    const int len_psi_grid = int(log2(N_timesteps)) + 1;;
    int recorder_psi_config[len_psi_grid];
    for (int ii=0; ii<len_psi_grid; ii++){recorder_psi_config[ii] = 0;}
    int recorder_S_basin[len_psi_grid];
    for (int ii=0; ii<len_psi_grid; ii++){recorder_S_basin[ii] = 0;}
    int recorder_E_basin[len_psi_grid];
    for (int ii=0; ii<len_psi_grid; ii++){recorder_E_basin[ii] = 0;}
    double total_time_in_E_basin = 0.0;
    double total_time_in_S_basin = 0.0;

    // ========================================================================
    // Aging functions ========================================================
    // ========================================================================

    int S_basin_pointer_1 = 0;  // Points to the current location on grid 1
    int S_basin_pointer_2 = 0;  // Points to the current location on grid 2
    int S_basin_index = 1;      // Indexes the last basin the tracer was/is in

    int E_basin_pointer_1 = 0;
    int E_basin_pointer_2 = 0;
    int E_basin_index = 1;

    // Fill the time grids for pi
    double time_pi_grid_1[n_samp];
    fill_pi_grid_1(time_pi_grid_1, N_timesteps, DW, n_samp);
    double time_pi_grid_2[n_samp];
    fill_pi_grid_2(time_pi_grid_2, time_pi_grid_1, DW, n_samp);

    // Initialize the recorder grids for pi
    int recorder_pi_basin_E_1[n_samp];
    int recorder_pi_basin_E_2[n_samp];
    int recorder_pi_basin_S_1[n_samp];
    int recorder_pi_basin_S_2[n_samp];

    // ========================================================================
    // Save the grids =========================================================
    // ========================================================================

    // Special case for index == 0. This only occurs once. Here, we save the
    // grids themselves to disk.
    if (index == 0)
    {
        dump_grid_to_disk(file_dump_loc, "energy_grid", time_energy_grid,
            n_samp);
        dump_grid_to_disk(file_dump_loc, "Pi_basin_grid_1", time_pi_grid_1,
            n_samp);
        dump_grid_to_disk(file_dump_loc, "Pi_basin_grid_2", time_pi_grid_2,
            n_samp);
    }


    // ========================================================================
    // Run the simulation =====================================================
    // ========================================================================

    while (energy_pointer < n_samp)
    {
        // Step 1: get the neighboring energies by filling the relevant object
        get_neighboring_energies(config, energy_arr, neighboring_energies,
            N_spins);

        // Step 2: get the exit rates and dE values
        get_exit_rates(current_energy, beta, neighboring_energies, exit_rates, 
            delta_E, N_spins);
        total_exit_rate = 0.0;
        for (int ii=0; ii<N_spins; ii++){total_exit_rate += exit_rates[ii];}

        // Step 3: initialize an exponential distribution to sample the
        // waiting time from
        std::exponential_distribution<double> tmp_exp_dist(total_exit_rate);
        waiting_time = tmp_exp_dist(generator);

        // At this point: we understand that the tracer is in the configuration
        // it was in at the beginning of this step, and that it is in that
        // configuration for `waiting_time` time, and the new current_time
        // after stepping is `current_time` + `waiting_time`.

        // Step 4: update the current time of the simulation clock
        current_time += waiting_time;

        // Step 5: step to the next state and store the proposed (new) energy
        step_next_state_(config, exit_rates, total_exit_rate, N_spins,
            generator);
        proposed_energy = energy_arr[binary_vector_to_int(config, N_spins)];


        // --------------------------------------------------------------------
        // Append the energy recorder -----------------------------------------
        // --------------------------------------------------------------------

        append_energy_tracker_(current_energy, current_time, time_energy_grid,
            recorder_energy_grid, &energy_pointer, n_samp);


        // --------------------------------------------------------------------
        // Append the psi recorders -------------------------------------------
        // --------------------------------------------------------------------

        // By definition, the psi_config recorder should be appended at every
        // step in the Gillespie algorithm since the tracer changes
        // configurations every time, and has spent waiting_time time in the
        // configuration
        append_psi_recorder_(recorder_psi_config, waiting_time, len_psi_grid);

        // If we're currently in a basin (current energy is less than a
        // threshold), then the tracer has been in that basin for waiting_time
        // time during this step.
        if (current_energy < thresh_S){total_time_in_S_basin += waiting_time;}
        if (current_energy < thresh_E){total_time_in_E_basin += waiting_time;}

        // Now, we check if the proposed energy is greater than the threshold,
        // since in that case we would have to append the recorder and reset
        // the counter
        if (proposed_energy >= thresh_S)
        {
            append_psi_recorder_(recorder_S_basin, total_time_in_S_basin,
                len_psi_grid);
            total_time_in_S_basin = 0.0;
        }
        if (proposed_energy >= thresh_E)
        {
            append_psi_recorder_(recorder_E_basin, total_time_in_E_basin,
                len_psi_grid);
            total_time_in_E_basin = 0.0;
        }


        // --------------------------------------------------------------------
        // Append the aging function recorders --------------------------------
        // --------------------------------------------------------------------

        // If we just entered a basin, update the basin index
        if (current_energy > thresh_E and proposed_energy <= thresh_E)
        {
            E_basin_index += 1;
        }
        if (current_energy > thresh_S and proposed_energy <= thresh_S)
        {
            S_basin_index += 1;
        }

        E_basin_pointer_1 = append_Pi_grid_recorder_(current_time,
            current_energy, thresh_E, E_basin_index, E_basin_pointer_1, n_samp,
            time_pi_grid_1, recorder_pi_basin_E_1);
        E_basin_pointer_2 = append_Pi_grid_recorder_(current_time,
            current_energy, thresh_E, E_basin_index, E_basin_pointer_2, n_samp,
            time_pi_grid_2, recorder_pi_basin_E_2);
        S_basin_pointer_1 = append_Pi_grid_recorder_(current_time,
            current_energy, thresh_E, S_basin_index, S_basin_pointer_1, n_samp,
            time_pi_grid_1, recorder_pi_basin_S_1);
        S_basin_pointer_2 = append_Pi_grid_recorder_(current_time,
            current_energy, thresh_E, S_basin_index, S_basin_pointer_2, n_samp,
            time_pi_grid_2, recorder_pi_basin_S_2);

        current_energy = proposed_energy;
    }

    dump_result_to_disk(file_dump_loc, "energy", recorder_energy_grid, n_samp,
        index);
    dump_result_ints_to_disk(file_dump_loc, "psi_c", recorder_psi_config,
        len_psi_grid, index);
    dump_result_ints_to_disk(file_dump_loc, "psi_b_S", recorder_S_basin,
        len_psi_grid, index);
    dump_result_ints_to_disk(file_dump_loc, "psi_b_E", recorder_E_basin,
        len_psi_grid, index);
    dump_result_ints_to_disk(file_dump_loc, "pi_b_S_1", recorder_pi_basin_S_1,
        n_samp, index);
    dump_result_ints_to_disk(file_dump_loc, "pi_b_S_2", recorder_pi_basin_S_2,
        n_samp, index);
    dump_result_ints_to_disk(file_dump_loc, "pi_b_E_1", recorder_pi_basin_E_1,
        n_samp, index);
    dump_result_ints_to_disk(file_dump_loc, "pi_b_E_2", recorder_pi_basin_E_2,
        n_samp, index);

    delete[] energy_arr;
}
