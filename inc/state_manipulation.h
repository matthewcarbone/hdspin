// This library contains all of the utilities required to manipulate
// the state objects (integer representations of spin systems)
// using arbitrary precision libraries.

#include "ArbitraryPrecision/ap/ap.hpp"

#ifndef STATE_MANIPULATION_H
#define STATE_MANIPULATION_H

// Default value for the arbitrary precision is 128
// Meaning we can have up to 128 spins
#ifndef PRECISON
#define PRECISON 128
#endif

namespace state_manipulation
{
    void get_neighbors_(ap_uint<PRECISON> *, ap_uint<PRECISON>,
        int);
    void arbitrary_precision_integer_from_int_array_(
        const int *, const int, ap_uint<PRECISON> &);
    void int_array_from_arbitrary_precision_integer_(
        int *, const int, ap_uint<PRECISON>);
}


#endif
