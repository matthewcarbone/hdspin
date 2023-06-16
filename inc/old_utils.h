/* Standalone utilities and helper functions.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include <random>

/* Computes the integer power with provided base and exponent. Note that both
base and exponent must be non-negative long long values, as this uses bitwise
manipulation to perform the exponentiation. */
long long ipow(long long, long long);

/* Converts a configuration (array) of spins (binary numbers, 0 & 1) into an
integer representation by assuming the vector is the binary representation of
the integer. Uses bitwise manipulations. */
long long binary_vector_to_int(const int *, const int);

/* Prints the configuration followed by the energy of that configuration. */
void print_config_and_energy(const double *, const int, const double);

/* Finds the minimum element in the provided array. Returns the index of that
element. */
int min_element(const double *, const int);

/* In-place flips the spin of the provided array at the location specified. */
void _helper_flip_spin_(int *, const int);

void load_long_long_grid_(std::vector<long long> &, const std::string);

double iterative_mean(double, double, long long);

double iterative_S(double, double, double, double);

double var_from_S(double, long long);

#endif
