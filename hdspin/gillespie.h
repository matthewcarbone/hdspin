#ifndef GILLESPIE
#define GILLESPIE

#include <vector>
#include <map>
#include <random>

void gillespie(const std::string file_dump_loc, const long int N_timesteps,
    const int N_spins, const double beta, const double beta_critical,
    const int landscape);

#endif

