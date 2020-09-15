#include <iostream>

#ifndef STRUCTURE_UTILS_H
#define STRUCTURE_UTILS_H


struct FileNames
{
    std::string energy;
    std::string psi_config;
    std::string aging_config_1;
    std::string aging_config_2;
    std::string ii_str;
    std::string grids_directory;
};

struct RuntimeParameters
{
    int log_N_timesteps;
    int N_spins;
    double beta;
    double beta_critical;
    int landscape;
    double energetic_threshold;
    double entropic_attractor;
};

FileNames get_filenames(const int, const std::string, const std::string);

RuntimeParameters get_runtime_params(const int, const int, const double,
    const double, const int);

#endif
