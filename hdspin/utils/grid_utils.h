#include <map>
#include <vector>
#include <iostream>

#ifndef GRID_UTILS_H
#define GRID_UTILS_H

void fill_pyLogspace(double *arr, const double start, const double stop,
    const int num, const double base);

void fill_pi_grid_1(double *arr, const double nMC,
    const double dw, const int n_samp);

void fill_pi_grid_2(double *arr, const double *pi_grid_1, const double dw,
    const int n_samp);

#endif
