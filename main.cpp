#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <chrono>
#include <omp.h>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <chrono>
#include <cassert>

// #include "cpprem/utils/filesystem_utils.h"
// #include "cpprem/utils/grid_utils.h"
// #include "cpprem/utils/initialization_utils.h"
#include "hdspin/gillespie.h"

// Defaults
// const double BETA_CRITICAL_EREM = 1.0;
// const double BETA_CRITICAL_REM = sqrt(2.0 * log(2.0));  // 1.0 / ~0.849
const int N_SAMP = 100;  // X log-spaced points on each sampling grid
// const int PRINT_EVERY_STANDARD = 1e5;
const int PRINT_EVERY_GILLESPIE = 1e5;
const double DW = 0.5;


/**
 * Prints info about the arguments and if needed, throws errors if the args
 * are invalid.
 */
void print_arg_info(const int dynamics_switch, const std::string file_dump_loc,
    const int n_spins, const double beta, const double bc,
    const double log10nmc, const long int nmc, const int landscape_switch,
    const double thresh_S, const double thresh_E)
{

    // Dynamics ---------------------------------------------------------------
    if (dynamics_switch == 0){printf("\tDynamics is STANDARD\n");}
    else if (dynamics_switch == 1)
    {
        printf("\tDynamics is GILLESPIE REJECTIONLESS\n");
    }
    else
    {
        throw std::runtime_error("\tInvalid DYNAMICS_SWITCH (1st arg)\n");
    }

    // File location ----------------------------------------------------------
    printf("\tFile dump location is '%s'\n", file_dump_loc.c_str());

    // Spins ------------------------------------------------------------------
    assert(n_spins > 0);
    printf("\tNumber of spins is %i\n", n_spins);

    // Inverse temperature ----------------------------------------------------
    assert(beta > 0);
    assert(bc > 0);
    printf("\tBeta is %.02f and critical beta is %.02f\n", beta, bc);

    // Timesteps --------------------------------------------------------------
    printf("\tLog10 NMC is %.02f => total %li timesteps\n", log10nmc, nmc);

    // Energy landscape -------------------------------------------------------
    if (landscape_switch == 0){printf("\tEnergy landscape is EREM\n");}
    else if (landscape_switch == 1){printf("\tEnergy landscape is REM\n");}
    else
    {
        throw std::runtime_error("Invalid LANDSCAPE_SWITCH (6th arg)\n");
    }

    // Energy thresholds ------------------------------------------------------
    printf("\tEntropic threshold is %.02f\n", thresh_S);
    printf("\tEnergetic threshold is %.02f\n", thresh_E);
}


/**
 * Main method 
 * -----------
 * Arguments are as follows:
 *      1. Dynamics switch
 *          = 0 if standard dynamics
 *          = 1 if Gillespie dynamics
 *      2. File dump location
 *          Indexes the location of where all results will be saved. Usually
 *          should be something like "results/run_type"
 *      3. The number of spins.
 *      4. The inverse temperature.
 *      5. The critical inverse temperature.
 *      6. The log10 value of the number of timesteps. For instance, if the
 *         value is 2.45, then then number of timesteps is int(10^2.45). This
 *         serves the purpose of the max time in rejectionless (Gillespie)
 *         simulations.
 *      7. Landscape switch
 *          = 0 for the EREM model
 *          = 1 for the REM model
 *      8. The number of independent serial runs to compute. The general idea
 *         is that each instance of this code will run on a single CPU and can
 *         perform multiple simulations before terminating.
 *      9. The index of the CPU being used in this computaiton. Used only in
 *         the offset in saving file names.
 */
int main(int argc, char const *argv[])
{ 
    printf("Program name is: %s\n", argv[0]); 
    if(argc != 10)  // = 1 + number of arguments in doc string above 
    { 
        throw std::runtime_error("Invalid number of command line args\n");
    }

    // Initialize all the command line arguemnts ------------------------------
    const int DYNAMICS_SWITCH = atoi(argv[1]);
    const std::string FILE_DUMP_LOC = argv[2];
    const int N_SPINS = atoi(argv[3]);
    const double BETA = atof(argv[4]);
    const double BETA_CRITICAL = atof(argv[5]);
    const long double log10nMC = atof(argv[6]);
    const long double base = 10.0;
    const long double double_NMC = pow(base, log10nMC);
    const long int NMC = ((long int) double_NMC);
    const int LANDSCAPE_SWITCH = atoi(argv[7]);
    const int N_SERIAL_RUNS = atoi(argv[8]);
    const int CPU_INDEX = atoi(argv[9]);

    std::cout << NMC << std::endl;

    // const int N_CONFIGS = int(pow(2, N_SPINS));
    // double *energy_arr = new double[N_CONFIGS];
    // initialize_energy_mapping_exponential_arr(energy_arr, N_SPINS, 1.0);
    // dump_arr_of_doubles_to_disk("results/e.txt", energy_arr, N_CONFIGS);
    // delete[] energy_arr;

    // Initialize the thresholds
    double THRESH_S, THRESH_E;
    if (LANDSCAPE_SWITCH == 0)  // EREM
    {
        THRESH_S = log((2.0 * BETA_CRITICAL - BETA) / BETA_CRITICAL)
            / (BETA - BETA_CRITICAL);
        THRESH_E = -log(N_SPINS) / BETA_CRITICAL;
    }
    else  // REM
    {
        THRESH_S = -N_SPINS * BETA / 2.0;
        THRESH_E = -sqrt(2.0 * N_SPINS * log(N_SPINS));
    }

    print_arg_info(DYNAMICS_SWITCH, FILE_DUMP_LOC, N_SPINS, BETA,
        BETA_CRITICAL, log10nMC, NMC, LANDSCAPE_SWITCH, THRESH_S, THRESH_E);

    // Begin the simulations --------------------------------------------------
    int run_idx;
    std::vector<double> elapsed_times;
    for (int serial_index=0; serial_index<N_SERIAL_RUNS; serial_index++)
    {
        run_idx = CPU_INDEX * N_SERIAL_RUNS + serial_index;
        auto start = std::chrono::high_resolution_clock::now();

        if (DYNAMICS_SWITCH == 0)
        {
            // run_local_standard_dynamics(run_idx, FILE_DUMP_LOC, NMC, N_SPINS,
            //     N_SAMP, BETA, BETA_CRITICAL, LANDSCAPE_SWITCH, THRESH_S,
            //     THRESH_E, PRINT_EVERY_STANDARD, DW);
            throw std::runtime_error("Not implemented");
        } 
        else
        {
            gillespie(run_idx, FILE_DUMP_LOC, NMC, N_SPINS,
                N_SAMP, BETA, BETA_CRITICAL, LANDSCAPE_SWITCH, THRESH_S,
                THRESH_E, PRINT_EVERY_GILLESPIE, DW);
        }

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration_seconds = 
            std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        double duration_double_seconds = 
            std::chrono::duration<double>(duration_seconds).count();
        printf("idx %i DONE in %.02f s\n", run_idx, duration_double_seconds);
        fflush(stdout);
        elapsed_times.push_back(duration_double_seconds);
    }

    double mean = 0.0;
    for (int ii=0; ii<N_SERIAL_RUNS; ii++){mean += elapsed_times[ii];}
    mean = mean / N_SERIAL_RUNS;
    double var = 0.0;
    for (int ii=0; ii<N_SERIAL_RUNS; ii++)
    {
        var += pow((mean - elapsed_times[ii]), 2);
    }
    var = var / (N_SERIAL_RUNS - 1);
    const double std = sqrt(var);
    printf("\nbeta %.02f N %i : timing mean/spread: %.02f +/- %.02f\n",
        BETA, N_SPINS, mean, std);

    return 0; 
} 
