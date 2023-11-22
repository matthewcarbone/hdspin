#ifndef MAIN_UTILS_H
#define MAIN_UTILS_H

#include "utils.h"
#include "spin.h"
#include "obs1.h"
#include "obs2.h"


// TODO: add documentation

namespace main_utils
{


/**
 * @brief [brief description]
 * @details [long description]
 * @return [description]
 */
void update_parameters_(utils::SimulationParameters* p);


void save_and_log_config(const utils::SimulationParameters p);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param p [description]
 */
void initialize_grids_and_directories(const utils::SimulationParameters p);


/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param params [description]
 */
void auto_determine_dynamics_(utils::SimulationParameters* params);


/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param params [description]
 */
void execute_process_pool(const utils::SimulationParameters params);

}

#endif
