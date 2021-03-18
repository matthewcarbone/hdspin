/* Auxiliary files general utilities.
 *
 * Matthew Carbone, Columbia University 2020
 */

#include <assert.h>
#include <vector>
#include <math.h>
#include <random>
#include <iostream>     // std::cout
#include <cstring>
#include <fstream>

#include "Utils/utils.h"


long long ipow(long long base, long long exp)
{
    assert(base > 0);
    assert(exp > 0);
    long long result = 1;
    for (;;)
    {
        if (exp & 1){result *= base;}
        exp >>= 1;
        if (!exp){break;}
        base *= base;
    }
    return result;
}

long long binary_vector_to_int(const int *config, const int N)
{
    long long res = 0;
    for (int ii=0; ii<N; ii++){res = res << 1LL | config[ii];} 
    return res;
}

void print_config_and_energy(const double *config, const int N,
    const double energy)
{
    for (int ii=0; ii<N; ii++)
    {
        std::cout << config[ii] << " ";
    }
    std::cout << energy << std::endl;
}


int min_element(const double *arr, const int N)
{
    int min_el = 0;
    for (int ii=1; ii<N; ii++)
    {
        if (arr[ii] < arr[min_el]){min_el = ii;}
    }
    return min_el;
}


void _helper_flip_spin_(int *cfg, const int idx)
{
    if (cfg[idx] == 0){cfg[idx] = 1;}
    else{cfg[idx] = 0;}
}


void _helper_calculate_neighboring_energies_(int *cfg,
    const double *energy_arr, const int N, double *neighboring_energies)
{
    for (int ii=0; ii<N; ii++)
    {
        // Flip the ii'th spin
        _helper_flip_spin_(cfg, ii);

        // Collect the energy of the spin_config
        neighboring_energies[ii] =
            energy_arr[binary_vector_to_int(cfg, N)];

        // Flip the ii'th spin back
        _helper_flip_spin_(cfg, ii);
    }
}

void load_long_long_grid_(std::vector<long long> &grid, const std::string loc)
{
    std::ifstream myfile (loc);
    std::string line;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            grid.push_back(stoll(line));
        }
        myfile.close();
    }
}

long double iterative_mean(long double mu0, long double x1, long long n)
{
    return mu0 + (x1 - mu0) / ((long double) n);
}


// S_n is a helper function for calculating the variance, such that
// var = S_n / n.
long double iterative_S(long double mu0, long double mu1, long double x1,
    long double S0)
{
    return S0 + (x1 - mu0) * (x1 - mu1);
}

long double var_from_S(long double S, long long n)
{
    assert(n >= 0);
    if (n == 0)
    {
        return 0.0;
    }
    return S / ((long double) n);
}
