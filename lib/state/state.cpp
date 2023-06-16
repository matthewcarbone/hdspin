#include <vector>
#include "state_manipulation.h"
#include "ArbitraryPrecision/ap/ap.hpp"

ap_uint<PRECISON> _arbitrary_precision_integer_pow(
    const int base, const int exponent)
{
    if (exponent == 0)
    {
        ap_uint<PRECISON> val = 1;
        return val;   
    }
    ap_uint<PRECISON> val = base;
    for (int ii=0; ii<exponent - 1; ii++)
    {
        val *= base;
    }
    return val;
}

namespace state_manipulation
{
    void get_neighbors_(ap_uint<PRECISON> *neighbors, ap_uint<PRECISON> n,
        int bitLength)
    {   
        for (int b = 0; b < bitLength; b++)
        {
            ap_uint<PRECISON> one = 1;
            neighbors[b] = n ^ (one << b);
        }
    }

    void arbitrary_precision_integer_from_int_array_(
        const int *config, const int N, ap_uint<PRECISON> &res)
    {
        res = 0;
        for (int ii=N-1; ii>=0; ii--)
        {
            const ap_uint<PRECISON> config_val = config[N - ii - 1];
            const ap_uint<PRECISON> power_val = _arbitrary_precision_integer_pow(2, ii);
            res += config_val * power_val;
        }
    }

    void int_array_from_arbitrary_precision_integer_(
        int *config, const int N, ap_uint<PRECISON> integer)
    {
        for (int ii=0; ii<N; ii++)
        {
            ap_uint<PRECISON> remainder = integer % 2;
            integer = integer / 2;
            config[N - ii - 1] = int(remainder);
        }
    }
}

