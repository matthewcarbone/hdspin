#include <iostream>

#ifndef STRUCTURE_UTILS_H
#define STRUCTURE_UTILS_H


struct Vals
{
    long long int_rep, int_rep_IS;
    long double energy, energy_IS;
};

struct FileNames
{
    std::string energy;
    std::string psi_config;
    std::string psi_basin;
    std::string aging_config_1;
    std::string aging_config_2;
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
    int dynamics;
    double energetic_threshold;
    double entropic_attractor;
};


RuntimeParameters get_runtime_params(const int, const int, const double,
    const double, const int);


#endif
