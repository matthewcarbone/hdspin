#include <iostream>
#include <vector>

// #include "ArbitraryPrecision/ap/ap.hpp"
#include "state_manipulation.h"

int main(int argc, char const *argv[])
{
    ap_uint<10000> inp = 0; // copy-initialization.
    std::vector<ap_uint<10000>> x = get_neighbors(inp, 32);

    for (int ii=0; ii<x.size(); ii++)
    {
        std::cout << x[ii] << std::endl;
    }
    return 0;
}
