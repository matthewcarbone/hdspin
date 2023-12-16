#ifndef PROCESSING_UTILS_H
#define PROCESSING_UTILS_H

#include "utils.h"

// std::vector<json> read_all_results();

// /**
//  * @brief Calculates all of the statistics of a standard result
//  * @details A standard result is given by one in which each tracer has a single
//  * vector of some length. The results that are calculated are the mean, median,
//  * standard deviation and standard error.
//  * 
//  * @param results A vector of loaded json objects
//  * @param key The key of the objects to access, e.g. "energy"
//  * @return A json object containing the results
//  */
// json get_standard_statistics(const std::vector<json> results, const std::string key);


namespace processing_utils 
{

/**
 * @brief [brief description]
 * @details [long description]
 */
void postprocess();

}

#endif
