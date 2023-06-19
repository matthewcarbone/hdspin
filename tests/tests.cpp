#include "catch2/catch_test_macros.hpp"  // New library version of Catch2
// #include "inc/Catch2/catch.hpp"  // Old header-only version of Catch2
#include "ArbitraryPrecision/ap/ap.hpp"

// Import the tests themselves
#include "test_utils.h"
#include "test_energy_mapping.h"


TEST_CASE("Test arbitrary precision interconversion", "[arbitrary_precision]")
{
    for (unsigned int seed=1; seed<11; seed++)
    {
        for (unsigned int arr_size=1; arr_size<11; arr_size++)
        {
            const unsigned int actual_arr_size = arr_size * 12;
            const unsigned int actual_seed = seed * 101;
            REQUIRE(test_array_to_int_to_array_conversion(actual_seed, actual_arr_size));
        }
    }
}

TEST_CASE("Test massive arbitrary precision (arr_size==PRECISON) nearest neighbors", "[arbitrary_precision]")
{
    const unsigned int seed = 1235;
    const unsigned int arr_size = PRECISON;

    // Note that the array size represents a number far larger than the
    // largest default precision long long integer
    REQUIRE(test_neighbors_correct(seed, arr_size));
}

TEST_CASE("Test energy mapping EREM sampling", "[energy_mapping]")
{
    for (int ii=1; ii<11; ii++)
    {
        const double beta_critical = (1.0 * ii - 0.5) * 0.5;
        REQUIRE(_test_energy_mapping_sampling_EREM_given_beta_critical(beta_critical));
    }
}

TEST_CASE("Test energy mapping REM sampling", "[energy_mapping]")
{
    for (int ii=1; ii<10; ii++)
    {
        REQUIRE(_test_energy_mapping_sampling_REM_given_N_spins(ii*10));
    }
}

TEST_CASE("Test small cache", "[energy_mapping]")
{
    for (int ii=1; ii<10; ii++)
    {
        REQUIRE(_test_small_cache(ii*10));
    } 
}

// int main(int argc, char const *argv[])
// {

//     return 0;
// }

// int main(int argc, char const *argv[])
// {

//     // const ap_uint<PRECISON> inp = *argv[1] - '0';
//     // const int nspins = *argv[2] - '0';

//     const int foo[96] = {
//         1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,1
//     };

//     ap_uint<PRECISON> res;
//     state_manipulation::arbitrary_precision_integer_from_int_array_(foo, 96, res);

//     int foo2[96];
//     state_manipulation::int_array_from_arbitrary_precision_integer_(foo2, 96, res);

//     for (int ii=0; ii<96; ii++)
//     {
//         std::cout << foo[ii] << " " << foo2[ii] << std::endl;
//     }


//     std::cout << res << std::endl;
// }
