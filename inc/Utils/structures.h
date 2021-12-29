#include <iostream>

#ifndef STRUCTURE_UTILS_H
#define STRUCTURE_UTILS_H


struct Vals
{
    long long int_rep;
    long long int_rep_IS;
    double energy;
    double energy_IS;
};

struct FileNames
{
    // Energy
    std::string energy, energy_IS, energy_avg_neighbors;
    std::string psi_config, psi_config_IS;
    std::string psi_basin_E, psi_basin_S, psi_basin_E_IS, psi_basin_S_IS;
    std::string aging_config, aging_config_IS;
    std::string aging_basin_E, aging_basin_E_IS, aging_basin_S, aging_basin_S_IS;
    std::string ridge_E_all, ridge_S_all, ridge_E_IS_all, ridge_S_IS_all;
    std::string unique_configs;
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
    std::string landscape;
    std::string dynamics;
    int divN;
    double energetic_threshold;
    double entropic_attractor;
    int memory;
    int memoryless_retain_last_energy;
    int max_ridges;
    bool valid_entropic_attractor;
    int n_tracers;
};

#endif
