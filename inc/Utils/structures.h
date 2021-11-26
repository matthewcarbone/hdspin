#include <iostream>

#ifndef STRUCTURE_UTILS_H
#define STRUCTURE_UTILS_H


struct Vals
{
    long long int_rep;
    long long int_rep_IS = -1;
    double energy;
    double energy_IS = 0.0;
};

struct FileNames
{
    // Energy
    std::string energy;
    std::string psi_config, psi_config_IS;
    std::string psi_basin_E, psi_basin_S, psi_basin_E_IS, psi_basin_S_IS;
    std::string aging_config_1;
    std::string aging_config_2;
    std::string aging_basin_1;
    std::string aging_basin_2;
    std::string ridge_E, ridge_S, ridge_E_IS, ridge_S_IS;
    std::string ridge_E_proxy_IS, ridge_S_proxy_IS;
    std::string ridge_E_all, ridge_S_all;
    std::string ii_str;
    std::string grids_directory;
};

struct RuntimeParameters
{
    int log_N_timesteps;
    int N_spins;
    long long N_configs;
    long long N_timesteps;
    double beta;
    double beta_critical;
    int landscape;
    int dynamics_flag;
    double energetic_threshold;
    double entropic_attractor;
    int memoryless;
    int max_ridges;
};

#endif
