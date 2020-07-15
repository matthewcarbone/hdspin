#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream


std::string filename(const std::string directory, const int identifier);

void dump_arr_of_doubles_to_disk(const std::string path,
    const double *v, const int N);

void dump_result_to_disk(const std::string base_loc,
    const std::string result_name, const double *v, const int N,
    const int index);

void dump_result_ints_to_disk(const std::string base_loc,
    const std::string result_name, const int *v, const int N,
    const int index);

void dump_grid_to_disk(const std::string base_loc,
    const std::string grid_name, const double *v, const int N);
