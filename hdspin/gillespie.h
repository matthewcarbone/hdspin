#ifndef GILLESPIE
#define GILLESPIE

#include <vector>
#include <map>
#include <random>

void gillespie(
    const int index,
    const std::string file_dump_loc,
    const long int N_timesteps,
    const int N_spins,
    const int n_samp,
    const double beta,
    const double beta_critical,
    const int landscape,
    const double thresh_S,
    const double thresh_E,
    const int print_every,
    const double DW);

#endif

