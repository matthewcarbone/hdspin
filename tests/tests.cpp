#include "catch2/catch_test_macros.hpp"  // New library version of Catch2
// #include "inc/Catch2/catch.hpp"  // Old header-only version of Catch2
#include "ArbitraryPrecision/ap/ap.hpp"

// Import the tests themselves
#include "test_utils.h"
#include "test_energy_mapping.h"
#include "test_spin.h"
#include "test_obs1.h"


TEST_CASE("Test arbitrary precision interconversion", "[arbitrary_precision]")
{
    std::cout << "PRECISON==" << PRECISON << std::endl;
    for (unsigned int seed=1; seed<11; seed++)
    {
        for (unsigned int arr_size=1; arr_size<11; arr_size++)
        {
            const unsigned int actual_arr_size = arr_size * 12;
            const unsigned int actual_seed = seed * 101;
            REQUIRE(test_utils::test_array_to_int_to_array_conversion(actual_seed, actual_arr_size));
        }
    }
}

TEST_CASE("Test massive arbitrary precision (arr_size==PRECISON) nearest neighbors", "[arbitrary_precision]")
{
    const unsigned int seed = 1235;

    // Note that the array size represents a number far larger than the
    // largest default precision long long integer
    REQUIRE(test_utils::test_neighbors_correct(seed, PRECISON));
}

TEST_CASE("Test flip_bit huge", "[arbitrary_precision]")
{
    REQUIRE(test_utils::test_flip_bit_small());
    REQUIRE(test_utils::test_flip_bit_huge(123, 5));
    REQUIRE(test_utils::test_flip_bit_huge(1234, 50));
    REQUIRE(test_utils::test_flip_bit_huge(12345, 100));
    REQUIRE(test_utils::test_flip_bit_huge(123456, PRECISON));
    REQUIRE(test_utils::test_flip_bit_big_number(PRECISON));
    REQUIRE(test_utils::test_flip_bit_big_number_self_consistent(PRECISON));
}

TEST_CASE("Test energy mapping EREM sampling", "[energy_mapping]")
{
    for (int ii=1; ii<11; ii++)
    {
        const double beta_critical = (1.0 * ii - 0.5) * 0.5;
        REQUIRE(test_energy_mapping::test_energy_mapping_sampling_EREM_given_beta_critical(beta_critical));
    }
}

TEST_CASE("Test energy mapping REM sampling", "[energy_mapping]")
{
    for (int ii=1; ii<10; ii++)
    {
        REQUIRE(test_energy_mapping::test_energy_mapping_sampling_REM_given_N_spins(ii*10));
    }
}

TEST_CASE("Test small cache", "[energy_mapping]")
{
    for (int ii=1; ii<10; ii++)
    {
        REQUIRE(test_energy_mapping::test_small_cache(ii*10));
    } 
}

TEST_CASE("Test memory -1", "[energy_mapping]")
{
    for (int ii=2; ii<12; ii++)
    {
        REQUIRE(test_energy_mapping::test_memory_minus_one(ii));
    }
}

TEST_CASE("Test massive AP LRU", "[energy_mapping]")
{
    int N_spins = 10;
    while (N_spins < PRECISON)
    {
        REQUIRE(test_energy_mapping::test_massive_AP_LRU(N_spins));
        N_spins += 10;
    }
}

TEST_CASE("Test spin emap", "[spin]")
{
    REQUIRE(test_spin::test_basic_lru_cache_with_spin());
}

TEST_CASE("Test inherent structure", "[spin]")
{
    // REQUIRE(test_inherent_structure());
    REQUIRE(test_spin::test_inherent_structure_min_is_min());
}

TEST_CASE("Test streaming median", "[obs1]")
{
    REQUIRE(test_obs1::test_streaming_median());
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
