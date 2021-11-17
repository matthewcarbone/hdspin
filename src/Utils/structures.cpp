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
    fnames.psi_config_IS = target_dir + "/" + ii_str + "_psi_config_IS.txt";

    fnames.psi_basin_E = target_dir + "/" + ii_str + "_psi_basin_E.txt";
    fnames.psi_basin_S = target_dir + "/" + ii_str + "_psi_basin_S.txt";
    fnames.psi_basin_E_IS = target_dir + "/" + ii_str + "_psi_basin_E_IS.txt";
    fnames.psi_basin_S_IS = target_dir + "/" + ii_str + "_psi_basin_S_IS.txt";

    fnames.aging_config_1 = target_dir + "/" + ii_str + "_pi1_config.txt";
    fnames.aging_config_2 = target_dir + "/" + ii_str + "_pi2_config.txt";

    fnames.aging_basin_1 = target_dir + "/" + ii_str + "_pi1_basin.txt";
    fnames.aging_basin_2 = target_dir + "/" + ii_str + "_pi2_basin.txt";

    fnames.ridge_E = target_dir + "/" + ii_str + "_ridge_E.txt";
    fnames.ridge_S = target_dir + "/" + ii_str + "_ridge_S.txt";
    fnames.ridge_E_IS = target_dir + "/" + ii_str + "_ridge_E_IS.txt";
    fnames.ridge_S_IS = target_dir + "/" + ii_str + "_ridge_S_IS.txt";

    fnames.ridge_E_proxy_IS = target_dir + "/" + ii_str + "_ridge_E_proxy_IS.txt";
    fnames.ridge_S_proxy_IS = target_dir + "/" + ii_str + "_ridge_S_proxy_IS.txt";

    fnames.ii_str = ii_str;
    fnames.grids_directory = grids_directory;
    return fnames;
}



RuntimeParameters get_runtime_params(const int log_N_timesteps,
    const int N_spins, const double beta, const double beta_critical,
    const int landscape, const int loop_dynamics, const int memoryless)
{
    double et, ea;
    RuntimeParameters params;
    params.log_N_timesteps = log_N_timesteps;
    params.N_spins = N_spins;
    params.beta = beta;
    params.beta_critical = beta_critical;
    params.landscape = landscape;
    params.loop_dynamics = loop_dynamics;
    params.memoryless = memoryless;

    if (landscape == 0) // EREM
    {
        et = -1.0 / beta_critical * log(N_spins);
        ea = 1.0 / (beta - beta_critical)
            * log((2.0 * beta_critical - beta) / beta_critical);

        if (beta >= 2.0 * beta_critical | beta <= beta_critical)
        {
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

    return params;
}
