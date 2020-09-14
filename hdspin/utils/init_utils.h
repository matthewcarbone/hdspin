#ifndef INIT_UTILS_H
#define INIT_UTILS_H

void initialize_spin_system(int *arr, const int N);

void initialize_energy_mapping_exponential_arr(double *,
    const unsigned long long, const double);

void initialize_energy_mapping_gaussian_arr(double *arr,
    const unsigned long long, const int N, const double bc);

#endif
