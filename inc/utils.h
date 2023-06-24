/**
 * This library contains all utilities required to manipulate the states in
 * the code, as well as general utils.
 */

#include <fstream>      // std::ofstream

#include "ArbitraryPrecision/ap/ap.hpp"
#include "Json/json.hpp"
using json = nlohmann::json;

#ifndef UTILS_H
#define UTILS_H

// Default value for the arbitrary precision is 128
// Meaning we can have up to 128 spins
// This must be defined at compile time using e.g. -DPRECISON=12345
// As per the Arbitrary Precision docs, this should probably be a power of 2!
#ifndef PRECISON
#define PRECISON 128
#endif

// The smoke test is a faster version of the tests
// this modifies a few of the tests that otherwise might take a while
// Set -DSMOKE=OFF when using GH actions to run tests.
// Otherwise the default is ON
#ifndef SMOKE
#define SMOKE 1
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

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param long [description]
 * @param long [description]
 * 
 * @return [description]
 */
long long ipow(long long base, long long exp);


// These functions are all defined in lib/state/state.cpp
namespace state
{

    /**
     * @brief Gets the neighboring states given a representation
     * @details Using binary operations (bit shifts) and arbitrary precision
     * integers, finds all neighbors in the binary bit-shift space.
     * 
     * @param n ap_uint<PRECISON> * Pointer to an array of neighbors to be populated
     * @param bitLength ap_uint<PRECISON> The binary representation of the state of which we
     * want to find the neighbors of
     * @param int The number of spins (bitlength)
     */
    void get_neighbors_(ap_uint<PRECISON> *neighbors, ap_uint<PRECISON> n, int bitLength);

    /**
     * @brief Flips a specific spin in the bit representation
     * 
     * @param state ap_uint<PRECISON> The current state
     * @param k int The bit to flip
     * 
     * @return ap_uint<PRECISON>
     */
    ap_uint<PRECISON> flip_bit(const ap_uint<PRECISON> state, const int k, const int bitLength);

    /**
     * @brief Converts an integer array to an arbitrary precision integer
     * @details Using the Arbitrary Precision library, converts a binary integer
     * array to an arbitrary precision integer
     * 
     * @param const int * Pointer to the binary array
     * @param const int The total number of spins
     * @param ap_uint<PRECISON> & Memory address of the integer to fill
     */
    void arbitrary_precision_integer_from_int_array_(const int *, const int, ap_uint<PRECISON> &);

    /**
     * @brief Reverses that of arbitrary_precision_integer_from_int_array_
     * @details Using the Arbitrary Precision library, converts an an arbitrary
     * precision integer to a binary integer array
     * 
     * @param int * Pointer to the binary array to fill
     * @param const int The number of spins
     * @param const ap_uint<PRECISON> The integer to compute
     */
    void int_array_from_arbitrary_precision_integer_(int *, const int, const ap_uint<PRECISON> &);

}

struct Counter
{
    unsigned long long rejections = 0;
    unsigned long long acceptances = 0;
};

namespace parameters
{

    // The properties of a state. Returned as a function of the states's
    // integer representation
    struct StateProperties
    {
        ap_uint<PRECISON> state;
        double energy;
    };

    struct FileNames
    {
        // Energy
        std::string energy, energy_IS, energy_avg_neighbors;
        std::string psi_config, psi_config_IS;
        std::string psi_basin_E, psi_basin_S, psi_basin_E_IS, psi_basin_S_IS;
        std::string aging_config, aging_config_IS;
        std::string aging_basin_E, aging_basin_E_IS, aging_basin_S, aging_basin_S_IS;
        std::string ridge_E_all, ridge_S_all, ridge_E_IS_all, ridge_S_IS_all;
        std::string unique_configs;
        std::string ii_str;
        std::string grids_directory;
    };

    struct SimulationParameters
    {
        int log10_N_timesteps;
        long long N_timesteps;
        unsigned int N_spins;
        std::string landscape;
        long long memory;
        double beta;
        double beta_critical;    
        std::string dynamics;
        int n_tracers_per_MPI_rank;

        unsigned int seed;
        bool use_manual_seed = false;
        bool valid_entropic_attractor = true;

        double energetic_threshold;
        double entropic_attractor;

        int grid_size = 100;
        double dw = 0.5;
    };

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param inp [description]
     */
    void log_json(const json inp);

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param p [description]
     */
    void log_parameters(const SimulationParameters p);

    /**
     * @brief [brief description]
     * @details [long description]
     * @return [description]
     */
    SimulationParameters get_parameters(const json inp);

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param ii [description]
     * @return [description]
     */
    FileNames get_filenames(const int ii);

}

/**
 * @brief [brief description]
 * @details [long description]
 */
void make_directories();

namespace grids
{

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param log10_timesteps [description]
     * @param n_gridpoints [description]
     */
    void make_energy_grid_logspace(const int log10_timesteps, const int n_gridpoints);

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param log10_timesteps [description]
     * @param dw [description]
     * @param n_gridpoints [description]
     */
    void make_pi_grids(const int log10_timesteps, const double dw, const int n_gridpoints);

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param long [description]
     * @param loc [description]
     * @param grid_size [description]
     */
    void load_long_long_grid_(std::vector<long long> &grid, const std::string loc);

}


namespace time_utils
{
    /**
     * @brief [brief description]
     * @details [long description]
     * @return [description]
     */
    std::string get_datetime();

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param start [description]
     * @return [description]
     */
    int get_time_delta(const std::chrono::time_point<std::chrono::high_resolution_clock> start);
}

#endif
