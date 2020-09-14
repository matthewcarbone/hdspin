#ifndef GILLESPIE
#define GILLESPIE

#include <iostream>

#include "utils/grid_utils.h"

void gillespie(EnergyGrid &, PsiConfigCounter &, const int, const int,
    const double, const double, const int);

#endif

