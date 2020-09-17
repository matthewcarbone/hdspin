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

struct SystemInformation
{
    long long x, x_prev;
    long double e, e_prev;

    // The amount of time the system has spent in some configuration
    long double waiting_time = 0.0;
};

FileNames get_filenames(const int, const std::string, const std::string);

RuntimeParameters get_runtime_params(const int, const int, const double,
    const double, const int);

#endif
