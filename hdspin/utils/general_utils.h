#include <random>

#ifndef GENERAL_UTILS_H
#define GENERAL_UTILS_H

long long ipow(long long base, long long exp);

long long binary_vector_to_int(const int *config, const int N);

void flip_spin_(int *config, const int idx);

void get_neighboring_energies(
    int *config, const double *energy_arr, double *energies,
    const int N);

void get_exit_rates(const double current_energy, const double beta,
    const double *neighboring_energies, double *exit_rates, 
    double *delta_E, const int N);

void step_next_state_(int *config, const double *exit_rates,
    const double total_exit_rate, const int N, std::mt19937 generator);

void print_config_and_energy(const double *config, const int N,
    const double energy);

long long compute_inherent_structure(const int *, const double *,
    const int);

long long query_inherent_structure(const int, const int *,
    const double *, long long *);

#endif
