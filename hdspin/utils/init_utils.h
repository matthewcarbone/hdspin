
#ifndef INIT_UTILS_H
#define INIT_UTILS_H

void initialize_spin_system(int *arr, const int N);
// std::map<int, double> energy_mapping_exponential(const int N, const double bc);
void initialize_energy_mapping_exponential_arr(double *arr,
    const int N, const double bc);

void initialize_energy_mapping_gaussian_arr(double *arr,
    const int N, const double bc);

#endif
