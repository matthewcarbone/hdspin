#ifndef GILLESPIE
#define GILLESPIE

#include <iostream>

void gillespie(const std::string file_dump_loc,
    const std::string is_path_dump_loc, const long int N_timesteps,
    const int N_spins, const double beta, const double beta_critical,
    const int landscape);

#endif

