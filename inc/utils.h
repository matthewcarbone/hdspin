/**
 * This library contains all utilities required to manipulate the states in
 * the code, as well as general utils.
 */

#include "ArbitraryPrecision/ap/ap.hpp"

#ifndef UTILS_H
#define UTILS_H

// Default value for the arbitrary precision is 128
// Meaning we can have up to 128 spins
// This must be defined at compile time using e.g. -DPRECISON=12345
// As per the Arbitrary Precision docs, this should probably be a power of 2!
#ifndef PRECISON
#define PRECISON 128
#endif

/**
 * @brief Defines arbitrary precision integer powers
 * @details Using the Arbitrary Precision library, defines custom code for
 * performing arbitrary precision power operations using only multiplication
 * in for loops
 * 
 * @param const int The base
 * @param const int The exponent
 * 
 * @return The integral result of the power operation.
 */
ap_uint<PRECISON> arbitrary_precision_integer_pow(const int, const int);


// These functions are all defined in lib/state/state.cpp
namespace state
{

/**
 * @brief Gets the neighboring states given a representation
 * @details Using binary operations (bit shifts) and arbitrary precision
 * integers, finds all neighbors in the binary bit-shift space.
 * 
 * @param ap_uint<PRECISON> * Pointer to an array of neighbors to be populated
 * @param ap_uint<PRECISON> The binary representation of the state of which we
 * want to find the neighbors of
 * @param int The number of spins (bitlength)
 */
void get_neighbors_(ap_uint<PRECISON> *, ap_uint<PRECISON>, int);

/**
 * @brief Converts an integer array to an arbitrary precision integer
 * @details Using the Arbitrary Precision library, converts a binary integer
 * array to an arbitrary precision integer
 * 
 * @param const int * Pointer to the binary array
 * @param const int The total number of spins
 * @param ap_uint<PRECISON> & Memory address of the integer to fill
 */
void arbitrary_precision_integer_from_int_array_(
    const int *, const int, ap_uint<PRECISON> &);

/**
 * @brief Reverses that of arbitrary_precision_integer_from_int_array_
 * @details Using the Arbitrary Precision library, converts an an arbitrary
 * precision integer to a binary integer array
 * 
 * @param int * Pointer to the binary array to fill
 * @param const int The number of spins
 * @param const ap_uint<PRECISON> The integer to compute
 */
void int_array_from_arbitrary_precision_integer_(
    int *, const int, const ap_uint<PRECISON> &);

}

#endif
