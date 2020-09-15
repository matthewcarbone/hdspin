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
    double energetic_threshold, entropic_attractor;
    RuntimeParameters params;
    params.log_N_timesteps = log_N_timesteps;
    params.N_spins = N_spins;
    params.beta = beta;
    params.beta_critical = beta_critical;
    params.landscape = landscape;

    if (landscape == 0) // EREM
    {
        energetic_threshold = -1.0 / beta_critical * log(N_spins);
        entropic_attractor = 1.0 / (beta - beta_critical)
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
        energetic_threshold = -sqrt(2.0 * N_spins * log(N_spins));
        entropic_attractor = -N_spins * beta / 2.0;
    }
    else
    {
        throw std::runtime_error("Unknown landscape flag");
    }

    params.energetic_threshold = energetic_threshold;
    params.entropic_attractor = entropic_attractor;

    // Arguments pertaining to the job itself
    printf("log_N_timesteps = %i\n", log_N_timesteps);
    printf("N_spins = %i\n", N_spins);
    printf("beta = %.02f\n", beta);
    printf("beta_critical = %.02f\n", beta_critical);
    printf("landscape = %i (0=EREM, 1=REM)\n", landscape);
    printf("energetic threshold = %.02f\n", energetic_threshold);
    printf("entropic attractor = %.02f\n", entropic_attractor);

    return params;
}
