#include <cassert>
#include <iostream>
#include <stdexcept>
#include <math.h>

#include "Utils/structures.h"
#include "Utils/utils.h"


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
    fnames.aging_basin_1 = target_dir + "/" + ii_str + "_pi1_basin.txt";
    fnames.aging_basin_2 = target_dir + "/" + ii_str + "_pi2_basin.txt";
    fnames.rolling = target_dir + "/" + ii_str + "_rolling.txt";
    fnames.ii_str = ii_str;
    fnames.grids_directory = grids_directory;
    return fnames;
}


RuntimeParameters get_runtime_params(const int log_N_timesteps,
    const int N_spins, const double beta, const double beta_critical,
    const int landscape, const int loop_dynamics)
{
    double et, ea;
    RuntimeParameters params;
    params.log_N_timesteps = log_N_timesteps;
    params.N_spins = N_spins;
    params.beta = beta;
    params.beta_critical = beta_critical;
    params.landscape = landscape;
    params.loop_dynamics = loop_dynamics;

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
            ea = 1e16;  // Set purposefully invalid value instead of nan or inf
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

    params.N_configs = ipow(2, params.N_spins);
    assert(params.N_configs > 0);  // Check for overflow
    params.N_timesteps = ipow(10, params.log_N_timesteps);
    assert(params.N_timesteps > 0);

    // Arguments pertaining to the job itself
    printf("log_N_timesteps = %i\n", params.log_N_timesteps);
    printf("N_timesteps = %lli\n", params.N_timesteps);
    printf("N_spins = %i\n", params.N_spins);
    printf("N_configs = %lli\n", params.N_configs);
    printf("beta = %.02f\n", params.beta);
    printf("beta_critical = %.02f\n", params.beta_critical);
    printf("landscape = %i (0=EREM, 1=REM)\n", params.landscape);
    printf("energetic threshold = %.02f\n", params.energetic_threshold);
    printf("entropic attractor = %.02f\n", params.entropic_attractor);

    return params;
}
