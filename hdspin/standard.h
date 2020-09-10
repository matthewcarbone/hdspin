#ifndef STANDARD
#define STANDARD

#include <iostream>

#include "utils/grid_utils.h"

void standard(EnergyGrid &energy_grid, const long long N_timesteps,
    const int N_spins, const double beta, const double beta_critical,
    const int landscape);

#endif

