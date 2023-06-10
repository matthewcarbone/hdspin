// This library contains all of the utilities required to manipulate
// the state objects (integer representations of spin systems)
// using arbitrary precision libraries.

#include <vector>
#include "ArbitraryPrecision/ap/ap.hpp"
 

std::vector<ap_uint<10000>> get_neighbors(ap_uint<10000> n, int bitLength)
{
    std::vector<ap_uint<10000>> result;
    for (int b = 0; b < bitLength; b++) {
        result.push_back(n ^ (1 << b));
    }
    return result;
}
