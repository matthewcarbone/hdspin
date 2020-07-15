/* Auxiliary files for initializing the objects we need.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <string>
#include <random>

#include "init_utils.h"


/**
 * Initialize a binary random vector of N spins
 */
void initialize_spin_system(int *arr, const int N)
{
    // Seed the RNG
    unsigned int seed = std::random_device{}();
    std::mt19937 generator(seed);

    // Initialize the distribution
    std::bernoulli_distribution distribution;
    
    // Fill the array
    for (int ii=0; ii<N; ii++){arr[ii] = distribution(generator);}
}



/**
 * Initialize a mapping between the integer representation of a binary list
 * and the value of the energy it takes on for an exponential distribution
 * (EREM model). This previous version is a dictionary, but there's no reason
 * to do it when we can simply use constant lookup time instead of log(N)
std::map<int, double> energy_mapping_exponential(const int N, const double bc)
{
    std::map<int, double> m;

    // Seed the RNG
    unsigned int seed = std::random_device{}();
    std::mt19937 generator(seed);

    // C++ exponential random number generator has input of the form
    // lambda * e^{-lambda * x} for x > 0. For us, the scale parameter is bc.
    std::exponential_distribution<double> distribution(bc);

    // Initialize the exponential distribution
    for (int ii=0; ii<int(pow(2, N)); ii++)
    {
        // Note for the EREM model the energies are all negative, hence why
        // we take the negative sign here.
        m[ii] = -distribution(generator);
    }

    return m;
}
*/


/**
 * Fill an array with the exponentially-sampled values for the energy.
 */
void initialize_energy_mapping_exponential_arr(double *arr,
    const int N, const double bc)
{
    // Seed the RNG
    unsigned int seed = std::random_device{}();
    std::mt19937 generator(seed);

    // C++ exponential random number generator has input of the form
    // lambda * e^{-lambda * x} for x > 0. For us, the scale parameter is bc.
    std::exponential_distribution<double> distribution(bc);

    // Initialize the exponential distribution
    for (int ii=0; ii<int(pow(2, N)); ii++)
    {
        // Note for the EREM model the energies are all negative, hence why
        // we take the negative sign here.
        arr[ii] = -distribution(generator);
    }
}


/**
 * Fill an array with the exponentially-sampled values for the energy.
 */
void initialize_energy_mapping_gaussian_arr(double *arr,
    const int N, const double bc)
{
    // Seed the RNG
    unsigned int seed = std::random_device{}();
    std::mt19937 generator(seed);

    std::normal_distribution<double> distribution(0.0, sqrt(N));

    for (int ii=0; ii<int(pow(2, N)); ii++)
    {
        arr[ii] = distribution(generator);
    }
}

