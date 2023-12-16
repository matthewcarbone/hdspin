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

// This is the current commit hash hdspin is running on
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "UNKNOWN_COMMIT_HASH"
#endif

// Define the grids directories and some paths
// These are constants always at compile time and should never change
// Always relative to the working directory in which the code was run
#define DATA_PATH "data/"
#define GRID_PATH "grids/"
#define ENERGY_GRID_PATH "grids/energy.txt"
#define PI1_GRID_PATH "grids/pi1.txt"
#define PI2_GRID_PATH "grids/pi2.txt"
#define CONFIG_PATH "config.json"
#define RESULTS_PATH "results.json"
#define DIAGNOSTICS_PATH "diagnostics.json"

namespace utils
{

/**
* @brief Computes the average of a vector of doubles.
*/
double mean_vector(const std::vector<double> v);
double weighted_mean_vector(const std::vector<double> v, const std::vector<double> weights);

/**
* @brief Computes the variance of a vector of doubles.
*/
double variance_vector(const std::vector<double> v);

double median_vector(const std::vector<double> v);

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


/**
 * @brief Gets the neighboring states given a representation
 * @details Using binary operations (bit shifts) and arbitrary precision
 * integers, finds all neighbors in the binary bit-shift space.
 * 
 * @param neighbors ap_uint<PRECISON>* Pointer to an array of neighbors to be populated
 * @param n ap_uint<PRECISON> The binary representation of the state of which we
 * want to find the neighbors of
 * @param bitLength int
 */
void get_neighbors_(ap_uint<PRECISON>* neighbors, ap_uint<PRECISON> n, const unsigned int bitLength);

/**
 * @brief Flips a specific spin in the bit representation
 * 
 * @param state ap_uint<PRECISON> The current state
 * @param k int The bit to flip
 * 
 * @return ap_uint<PRECISON>
 */
ap_uint<PRECISON> flip_bit(const ap_uint<PRECISON> state, const unsigned int k, const unsigned int bitLength);

/**
 * @brief Converts an integer array to an arbitrary precision integer
 * @details Using the Arbitrary Precision library, converts a binary integer
 * array to an arbitrary precision integer
 * 
 * @param const int * Pointer to the binary array
 * @param const int The total number of spins
 * @param ap_uint<PRECISON> & Memory address of the integer to fill
 */
void arbitrary_precision_integer_from_int_array_(const unsigned int *, const unsigned int, ap_uint<PRECISON> &);

/**
 * @brief Reverses that of arbitrary_precision_integer_from_int_array_
 * @details Using the Arbitrary Precision library, converts an an arbitrary
 * precision integer to a binary integer array
 * 
 * @param int * Pointer to the binary array to fill
 * @param const int The number of spins
 * @param const ap_uint<PRECISON> The integer to compute
 */
void int_array_from_arbitrary_precision_integer_(unsigned int *, const unsigned int, const ap_uint<PRECISON> &);


/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param current_state [description]
 * @param int [description]
 * 
 * @return [description]
 */
std::string string_rep_from_arbitrary_precision_integer(const ap_uint<PRECISON> current_state, const unsigned int N);

// The properties of a state. Returned as a function of the states's
// integer representation
struct StateProperties
{
    ap_uint<PRECISON> state;
    double energy;
};

struct FileNames
{
    std::string ii_str, json_final;
};

struct SimulationParameters
{
    // None of these should have defaults
    unsigned int log10_N_timesteps;
    unsigned int N_spins;
    std::string landscape;
    double beta;

    // Some come along with defaults
    long long memory = pow(2, 25);
    std::string dynamics = "auto";
    unsigned int n_tracers = 100;
    unsigned int seed = 0;  // 0 is special, meaning no seed

    // Some defaults which are not required to be explicitly set by the user
    long long N_timesteps;
    double beta_critical;
    double energetic_threshold;
    double entropic_attractor;
    bool use_manual_seed = false;
    // bool valid_energetic_threshold = true;
    bool valid_entropic_attractor = true;
    int grid_size = 100;
    double dw = 0.5;
    bool calculate_inherent_structure_observables = false;
};

struct SimulationStatistics
{
    unsigned long long rejections = 0;
    unsigned long long acceptances = 0;
    unsigned long long total_steps = 0;
    double total_wall_time = 0.0;
    double total_waiting_time = 0.0;
};

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param inp [description]
 */
void print_json(const json jrep);

json simulation_parameters_to_json(const utils::SimulationParameters p);

void json_to_file(const json jrep, const std::string& filename);
void json_to_file_no_format(const json jrep, const std::string& filename);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param ii [description]
 * @return [description]
 */
FileNames get_filenames(const unsigned int ii);


void make_directories();

void cleanup_directories();

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
double get_time_delta(const std::chrono::time_point<std::chrono::high_resolution_clock> start);

}

#endif
