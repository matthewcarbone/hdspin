#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <random>

#include "utils.h"


void _fill_binary_vector(int *arr, const unsigned int arr_size,
    const unsigned int seed)
{
    std::mt19937 generator;
    // unsigned int seed = std::random_device{}();
    generator.seed(seed);
    std::bernoulli_distribution bernoulli_distribution;

    for (int ii=0; ii<arr_size; ii++)
    {
        arr[ii] = bernoulli_distribution(generator);
    }
}

/**
 * @brief Tests the conversion between an integer array and decimal
 * representation
 * @details First, fills an array of size arr_size with random 0's and 1's.
 * Next, converts that binary vector to decimal representation, then
 * converts it back. Then tests that the new vector has the same ordering
 * of 0's and 1's as the initial array.
 * 
 * @param int The random seed
 * @param int The size of the binary array
 * 
 * @return Returns true if all checks pass, false otherwise.
 */
bool test_array_to_int_to_array_conversion(
    const unsigned int seed, const unsigned int arr_size)
{
    int arr[arr_size];
    _fill_binary_vector(arr, arr_size, seed);

    ap_uint<PRECISON> res;
    state::arbitrary_precision_integer_from_int_array_(arr, arr_size, res);

    int arr2[arr_size];
    state::int_array_from_arbitrary_precision_integer_(arr2, arr_size, res);

    for (int ii=0; ii<arr_size; ii++)
    {
        if (arr[ii] != arr2[ii]){return false;}
    }
    return true;
}

/**
 * @brief Compares nearest neighbors computed form the vector and decimal
 * representations
 * 
 * @details This code executes the following steps:
 * 1. Creates a random binary vector based on the size of the array and
 * the seed.
 * 2. From that vector, flips each spin manually, computes the decimal
 * representation and saves that integer.
 * 3. From that vector, computes the decimal representation, then uses
 * bit-flip operations to find the neighbors in that representation.
 * 4. Compares the representations.
 * 
 * @param int The random seed
 * @param int The size of the binary array.
 * 
 * @return Returns true if all checks pass, false otherwise.
 */
bool test_neighbors_correct(const unsigned int seed, const unsigned int arr_size)
{
    int arr[arr_size];
    _fill_binary_vector(arr, arr_size, seed);

    // Get the neighbor values via flipping spins manually
    ap_uint<PRECISON> neighbors_method_1[arr_size];
    for (int ii=0; ii<arr_size; ii++)
    {
        // Flip the ii'th bit
        if (arr[ii] == 0){arr[ii] = 1;}
        else{arr[ii] = 0;}

        // Get the decimal representation of this binary vector
        state::arbitrary_precision_integer_from_int_array_(arr, arr_size, neighbors_method_1[ii]);

        // Flip the ii'th bit back
        if (arr[ii] == 0){arr[ii] = 1;}
        else{arr[ii] = 0;}
    }

    // Now we do the opposite: first convert the array to decimal
    // representation, and find its neighbors from that directly.
    ap_uint<PRECISON> decimal_rep;
    state::arbitrary_precision_integer_from_int_array_(arr, arr_size, decimal_rep);

    ap_uint<PRECISON> neighbors_method_2[arr_size];
    state::get_neighbors_(neighbors_method_2, decimal_rep, arr_size);

    // Compare!
    // Note they save the neighbors in opposite order!
    for (int ii=0; ii<arr_size; ii++)
    {
        const ap_uint<PRECISON> val1 = neighbors_method_1[ii];
        const ap_uint<PRECISON> val2 = neighbors_method_2[ii];
        // std::cout << val1 << " " << val2 << " " << std::string(val1) << " " << std::string(val2) << std::endl;
        if (val1 != val2){return false;}
    }

    return true;
}

bool test_flip_bit_small()
{
    ap_uint<PRECISON> thirteen = 13;
    ap_uint<PRECISON> v2 = state::flip_bit(thirteen, 1, 4);
    if (v2 == 9){return true;}
    return false;
}

bool test_flip_bit_huge(const unsigned int seed, const unsigned int arr_size)
{
    int arr[arr_size];
    _fill_binary_vector(arr, arr_size, seed);
    arr[0] = 1;
    arr[arr_size - 1] = 1;
    ap_uint<PRECISON> v1, v2, ap_rep;
    state::arbitrary_precision_integer_from_int_array_(arr, arr_size, ap_rep);

    for (int ii=0; ii<arr_size; ii++)
    {
        // Flip in the spin representation
        if (arr[ii] == 0){arr[ii] = 1;}
        else{arr[ii] = 0;}

        // Calculate the value
        state::arbitrary_precision_integer_from_int_array_(arr, arr_size, v1);

        // Flip back
        if (arr[ii] == 0){arr[ii] = 1;}
        else{arr[ii] = 0;}

        // Calculate the value by flipping bits
        v2 = state::flip_bit(ap_rep, ii, arr_size);

        // std::cout << v1 << " " << v2 << std::endl;
        if (v1 != v2){return false;}
    }

    return true;
}


bool test_flip_bit_big_number(const unsigned int arr_size)
{
    int arr[arr_size];
    for (int ii=0; ii<arr_size; ii++){arr[ii] = 1;}
    
    ap_uint<PRECISON> v1, v2, ap_rep, v2_prime;
    state::arbitrary_precision_integer_from_int_array_(arr, arr_size, ap_rep);

    for (int ii=0; ii<arr_size; ii++)
    {
        // Flip in the spin representation
        if (arr[ii] == 0){arr[ii] = 1;}
        else{arr[ii] = 0;}

        // Calculate the value
        state::arbitrary_precision_integer_from_int_array_(arr, arr_size, v1);

        // Flip back
        if (arr[ii] == 0){arr[ii] = 1;}
        else{arr[ii] = 0;}

        // Calculate the value by flipping bits
        v2 = state::flip_bit(ap_rep, ii, arr_size);

        // std::cout << v1 << " " << v2 << std::endl;
        if (v1 != v2){return false;}
    }

    return true;
}


bool test_flip_bit_big_number_self_consistent(const unsigned int arr_size)
{
    int arr[arr_size];
    for (int ii=0; ii<arr_size; ii++){arr[ii] = 1;}
    
    ap_uint<PRECISON> ap_rep, flipped, ap_rep_2;
    state::arbitrary_precision_integer_from_int_array_(arr, arr_size, ap_rep);

    for (int ii=0; ii<arr_size; ii++)
    {
        // Flip a bit
        flipped = state::flip_bit(ap_rep, ii, arr_size);

        // Flip it back
        ap_rep_2 = state::flip_bit(flipped, ii, arr_size);

        // std::cout << v1 << " " << v2 << std::endl;
        if (ap_rep != ap_rep_2){return false;}
    }

    return true;
}

#endif
