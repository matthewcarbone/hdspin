#include <iostream>
#include <stdexcept>
#include <math.h>

#include "structure_utils.h"


FileNames get_filenames(const int ii, const std::string target_dir,
    const std::string grids_directory)
{
    std::string ii_str = std::to_string(ii);
    ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');
    FileNames fnames;
    fnames.energy = target_dir + "/" + ii_str + "_energy.txt";
    fnames.psi_config = target_dir + "/" + ii_str + "_psi_config.txt";
    fnames.psi_basin = target_dir + "/" + ii_str + "_psi_basin.txt";
    fnames.aging_config_1 = target_dir + "/" + ii_str + "_pi1_config.txt";
    fnames.aging_config_2 = target_dir + "/" + ii_str + "_pi2_config.txt";
    fnames.ii_str = ii_str;
    fnames.grids_directory = grids_directory;
    return fnames;
}


RuntimeParameters get_runtime_params(const int log_N_timesteps,
    const int N_spins, const double beta, const double beta_critical,
    const int landscape)
{
    double et, ea;
    RuntimeParameters params;
    params.log_N_timesteps = log_N_timesteps;
    params.N_spins = N_spins;
    params.beta = beta;
    params.beta_critical = beta_critical;
    params.landscape = landscape;

    if (landscape == 0) // EREM
    {
        et = -1.0 / beta_critical * log(N_spins);
        ea = 1.0 / (beta - beta_critical)
            * log((2.0 * beta_critical - beta) / beta_critical);

        if (beta >= 2.0 * beta_critical | beta <= beta_critical)
        {
            printf(
                "---WARNING: beta restriction bc < b < 2bc not satisfied---\n"
            );
        }
    }
    else if (landscape == 1) // REM
    {
        et = -sqrt(2.0 * N_spins * log(N_spins));
        ea = -N_spins * beta / 2.0;
    }
    else
    {
        throw std::runtime_error("Unknown landscape flag");
    }

    params.energetic_threshold = et;
    params.entropic_attractor = ea;

    // Arguments pertaining to the job itself
    printf("log_N_timesteps = %i\n", log_N_timesteps);
    printf("N_spins = %i\n", N_spins);
    printf("beta = %.02f\n", beta);
    printf("beta_critical = %.02f\n", beta_critical);
    printf("landscape = %i (0=EREM, 1=REM)\n", landscape);
    printf("energetic threshold = %.02f\n", et);
    printf("entropic attractor = %.02f\n", ea);

    return params;
}


void update_basin_information(SystemInformation * system,
    const RuntimeParameters params, const double waiting_time)
{

    // We just entered a basin
    if (system->e < params.energetic_threshold &
        system->e_prev >= params.energetic_threshold)
    {
        system->basin_energy += 1;

        // This value is appended to the counter at every iteration of the
        // outer while loop. This is checked and reset by the basin tracker. If
        // its value is non-zero it will be appended to the counter. Else, it
        // will not be appended. This allows for compatibility with both
        // standard and Gillespie dynamics.
        system->tmp_t_basin_energy = system->t_basin_energy;

        system->t_basin_energy = 0.0;
    }

    // Regardless of if we exited the basin in the next step, we always
    // append the waiting time if the previous energy was below the threshold
    else if (system->e_prev < params.energetic_threshold)
    {
        system->t_basin_energy += waiting_time;
    }


    // Same update for the entroyp

    // We just entered a basin
    if (system->e < params.entropic_attractor &
        system->e_prev >= params.entropic_attractor)
    {
        system->basin_entropy += 1;
        system->tmp_t_basin_entropy = system->t_basin_entropy;
        system->t_basin_entropy = 0.0;
    }

    // Regardless of if we exited the basin in the next step, we always
    // append the waiting time if the previous energy was below the threshold
    else if (system->e_prev < params.entropic_attractor)
    {
        system->t_basin_entropy += waiting_time;
    }

}
