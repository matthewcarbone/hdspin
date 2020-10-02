#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "Utils/utils.h"

// TEST UTILS -----------------------------------------------------------------

TEST_CASE("Test ipow", "[utils]")
{
    REQUIRE(ipow(1, 7) == 1);
    REQUIRE(ipow(1, 40) == 1);
    REQUIRE(ipow(2, 4) == 16);
    REQUIRE(ipow(7, 21) == 558545864083284007);
    REQUIRE(ipow(7, 21) > 558545864083284006);
    REQUIRE(ipow(10, 18) == 1000000000000000000);
    REQUIRE(ipow(2, 36) == 68719476736);
}

TEST_CASE("Binary vector to integer conversion", "[utils]")
{
    const int arr1[18] = {0,1,0,0,1,0,1,0,0,1,0,0,1,1,1,0,0,1};
    REQUIRE(binary_vector_to_int(arr1, 18) == 76089);
    const int arr2[30] = 
        {1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,0,1,0,0,0,0,0,1,1,1,1};
    REQUIRE(binary_vector_to_int(arr2, 30) == 1057462799);
    const int arr3[10] = {0,0,0,0,0,0,0,0,0,0};
    REQUIRE(binary_vector_to_int(arr3, 10) == 0);
}

TEST_CASE("Minimum element", "[utils]")
{
    const double arr1[5] = {0.1234, 0.2, 0.3, 0.4, 0.5};
    REQUIRE(min_element(arr1, 5) == 0);
    const double arr2[5] = {0.1234, 0.2, 0.3, -0.4, 0.0};
    REQUIRE(min_element(arr2, 5) == 3);
    const double arr3[10] = {0.1234, 0.2, 0.3, -0.4, 0.0, 10, -10, 5, 6.5, 9};
    REQUIRE(min_element(arr3, 10) == 6);
}

TEST_CASE("Flip spin", "[utils]")
{
    int spins[8] = {0, 1, 1, 0, 1, 1, 0, 1};
    const int spins_copy[8] = {0, 1, 1, 0, 1, 1, 0, 1};

    for (int s=0; s<8; s++)
    {
        _helper_flip_spin_(spins, s);
        for (int ii=0; ii<8; ii++)
        {
            if (ii == s){REQUIRE(abs(spins[ii] - 1) == spins_copy[ii]);}
            else{REQUIRE(spins[ii] == spins_copy[ii]);}
        }
        _helper_flip_spin_(spins, s);   
    }
}

TEST_CASE("Neighboring energies", "[utils]")
{
    const double emap[4] = {0.123, -0.234, 0.345, -0.987};
    int config[2] = {0, 1};
    double ne[2];
    _helper_calculate_neighboring_energies_(config, emap, 2, ne);
    // config initial condition imples flipping:
    // spin 0 => [1, 1] => intrep = 3 => -0.987
    // spin 1 => [0, 0] => intrep = 0 => 0.123
    REQUIRE(ne[0] == -0.987);
    REQUIRE(ne[1] == 0.123);
}
