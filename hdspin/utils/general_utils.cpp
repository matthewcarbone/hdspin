/* Auxiliary files general utilities.
 *
 * Matthew Carbone, Columbia University 2020
 */

#include <vector>
#include <math.h>
#include <random>
#include <iostream>     // std::cout
#include <sstream>

#include "general_utils.h"

int binary_vector_to_int(const int *config, const int N)
{
    int res = 0;
    for (int ii=0; ii<N; ii++){res = res << 1 | config[ii];}    
    return res;
}


/**
 * Modifies the config vector in place by flipping the spin at the idx'th
 * entry.
 * @param idx int the index to flip the spin at.
 */
void flip_spin_(int *config, const int idx)
{
    if (config[idx] == 0){config[idx] = 1;}
    else{config[idx] = 0;}
}


/**
 * Get's the neighboring energies of the current configuration.
 */
void get_neighboring_energies(
    int *config, const double *energy_arr, double *energies,
    const int N)
{
    for (int ii=0; ii<N; ii++)
    {
        // Flip the ii'th spin
        flip_spin_(config, ii);

        // Collect the energy of the configuration
        energies[ii] = energy_arr[binary_vector_to_int(config, N)];

        // Flip the ii'th spin back
        flip_spin_(config, ii);
    }
}


/**
 * Get's the exit rates and delta_E.
 */
void get_exit_rates(const double current_energy, const double beta,
    const double *neighboring_energies, double *exit_rates, 
    double *delta_E, const int N)
{
    for (int ii=0; ii<N; ii++)
    {
        delta_E[ii] = neighboring_energies[ii] - current_energy;
        exit_rates[ii] = exp(-beta * delta_E[ii]);
        if (exit_rates[ii] > 1.0){exit_rates[ii] = 1.0;}
    }
    for (int ii=0; ii<N; ii++){exit_rates[ii] = exit_rates[ii] / N;}
}


void step_next_state_(int *config, const double *exit_rates,
    const double total_exit_rate, const int N, std::mt19937 generator)
{
    std::vector<double> normalized_exit_rates;
    for (int ii=0; ii<N; ii++)
    {
        normalized_exit_rates.push_back(exit_rates[ii] / total_exit_rate);
    }

    // Now, we make a choice of the spin to flip
    std::discrete_distribution<int> _dist(
        normalized_exit_rates.begin(), normalized_exit_rates.end());
    const int spin_to_flip = _dist(generator);
    flip_spin_(config, spin_to_flip);
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


/* Finds the lowest energy configuration via greedy search from a config
 * through its neighbors. Returns the config index of this inherent
 structure. */
int compute_inherent_structure(const int *config, const double *energy_arr,
    const int N)
{

    int min_el;

    // Make a copy of the config
    int config_copy[N];
    memcpy(config_copy, config, N*sizeof(int));

    double neighboring_energies[N];
    double current_energy;

    int config_int;

    // While we are not at the inherent structure, keep going
    while (true)
    {
        config_int = binary_vector_to_int(config_copy, N); 

        // Get the initial energies
        get_neighboring_energies(config_copy, energy_arr, neighboring_energies,
            N);

        // Get the current energy
        current_energy = energy_arr[config_int];

        // If we have not reached the lowest local energy which can be reached
        // by flipping a spin
        min_el = min_element(neighboring_energies, N);
        if (neighboring_energies[min_el] < current_energy)
        {
            flip_spin_(config_copy, min_el);
        }
        else
        {
            return binary_vector_to_int(config_copy, N);
        }
    }
}


std::string ZeroPadNumber(int num)
{
    std::stringstream ss;
    
    // the number is converted to string with the help of stringstream
    ss << num; 
    std::string ret;
    ss >> ret;
    
    // Append zero chars
    int str_length = ret.length();
    for (int i = 0; i < 7 - str_length; i++)
        ret = "0" + ret;
    return ret;
}
