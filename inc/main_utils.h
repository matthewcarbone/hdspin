#ifndef MAIN_UTILS_H
#define MAIN_UTILS_H

#include "utils.h"
#include "spin.h"
#include "obs1.h"
#include "obs2.h"



/**
 * The main_utils namespace provides an API for the entrypiont in main.cpp.
 * This way everything can be properly scoped and organized. main.cpp will only
 * call functions from this namespace.
 */
namespace main_utils
{


/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param p [description]
 */
void simulation_parameters_from_disk_(utils::SimulationParameters* p);


/**
 * @brief Completes the CLI-provided parameters
 * @details The user has full control over almost every part of the simulation,
 * but can optionally not provide explicit arguments for all of them. Additionally,
 * some parameters are explicit functions of what the user provides. This method
 * modifies the simulation parameters in-place.
 * 
 * @param p Simulation parameters
 */
void update_parameters_(utils::SimulationParameters* p);


/**
 * @brief Saves and logs the config file
 * @details Given a SimulationParameters struct, this first converts the
 * object to json format, then logs it (prints to std::cout), then saves it as
 * config.json in the working directory.
 * 
 * @param p Simulation parameters
 */
void save_and_log_config(const utils::SimulationParameters p);

/**
 * @brief Initializes directories and the grid utility files
 * @details hdspin requires an energy and "pi" grid. This function creates the
 * data and grids directories, then populates the grid directory with the
 * aforementioned files.
 * 
 * @param p Simulation parameters
 */
void initialize_grids_and_directories(const utils::SimulationParameters p);


/**
 * @brief Given "auto" as a dynamics type, automatically determines which
 * dynamics will be faster (gillespie or standard)
 * @details This is an MPI-aware function which uses all available tasks to get
 * some lightweight statistics on which simulation dynamics will be faster.
 * Specifically, it calculates the average walltime/simulation time (this is
 * necessary since gillespie can, by design, sometimes overshoot the maximum
 * simulation time). The faster simulator is then set in the simulation
 * parameters in-place.
 * 
 * @param p Simulation parameters
 */
void auto_determine_dynamics_(utils::SimulationParameters* params);


/**
 * @brief Executes the simulation
 * @details Utilizes a process pool: a single "master" controller dynamically
 * allocates tasks to workers until the total number of simulations is complete.
 * 
 * @param p Simulation parameters
 */
void execute_process_pool(const utils::SimulationParameters params);

}

#endif
