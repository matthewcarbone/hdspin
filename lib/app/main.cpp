#include <iostream>
#include <vector>

#include "../state_manipulation/state_manipulation.h"

// ap_uint<PRECISON> from_string(const char *string)
// {
//     *this = 0;
//     while( *string != '\0' ) {
//         *this *= 10;
//         *this += *string - '0';
//         string++;
//     }
// }

// int main(int argc, char const *argv[])
// {

//     // const ap_uint<PRECISON> inp = *argv[1] - '0';
//     // const int nspins = *argv[2] - '0';

//     const ap_uint<PRECISON> inp = atoi(argv[1]);
//     const int nspins = atoi(argv[2]);
//     std::cout << inp << " " << nspins << std::endl;

//     std::vector<ap_uint<PRECISON>> x = state_manipulation::get_neighbors(inp, nspins);

//     for (int ii=0; ii<x.size(); ii++)
//     {
//         std::cout << ii << " ~ " << x[ii] << std::endl;
//     }
//     return 0;
// }


int main(int argc, char const *argv[])
{

    // const ap_uint<PRECISON> inp = *argv[1] - '0';
    // const int nspins = *argv[2] - '0';

    const int foo[48] = {
        1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1
    };

    ap_uint<PRECISON> res;
    state_manipulation::arbitrary_precision_integer_from_int_array_(foo, 48, res);

    int foo2[48];
    state_manipulation::int_array_from_arbitrary_precision_integer_(foo2, 48, res);

    for (int ii=0; ii<48; ii++)
    {
        std::cout << foo[ii] << " " << foo2[ii] << std::endl;
    }


    std::cout << res << std::endl;
}
