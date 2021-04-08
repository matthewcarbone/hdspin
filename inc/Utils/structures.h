#include <iostream>

#ifndef STRUCTURE_UTILS_H
#define STRUCTURE_UTILS_H


struct Vals
{
    long long int_rep, int_rep_IS;
    double energy, energy_IS;
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
    std::string ii_str;
    std::string grids_directory;
};

FileNames get_filenames(const int, const std::string, const std::string);

struct RuntimeParameters
{
    int log_N_timesteps;
    long long N_timesteps;
    int N_spins;
    long long N_configs;
    double beta;
    double beta_critical;
    int landscape;
    double energetic_threshold;
    double entropic_attractor;
    int loop_dynamics;
};


RuntimeParameters get_runtime_params(const int, const int, const double,
    const double, const int, const int);


#endif
