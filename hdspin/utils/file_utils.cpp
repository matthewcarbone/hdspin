#include <string>       // std::string
#include <iostream>     // std::cout
#include <iomanip>      // std::setfill, std::setw
#include <sstream>      // std::stringstream
#include <vector>
#include <fstream>      // std::ofstream

#include "file_utils.h"

const int N_ZERO_PAD = 10;

/* Returns a string consisting of the combined filename. For example, for
 * directory = "results/energy" and identifier 10, it will return the string
 * "results/energy/000000010.txt", depending on the number of integers to
 * zero pad with.
 */
std::string filename(const std::string directory, const int identifier)
{
    std::ostringstream stringStream;
    stringStream << directory << "/"
        << std::setfill('0') << std::setw(N_ZERO_PAD) << identifier
        << ".txt";
    const std::string finalString = stringStream.str();
    return finalString;
}

/* Output a vector of doubles to disk, see
 * https://stackoverflow.com/questions/6406356
 * /how-to-write-vector-values-to-a-file
 */
void dump_arr_of_doubles_to_disk(const std::string path,
    const double *v, const int N)
{
    std::ofstream outFile(path);
    for (int ii=0; ii<N; ii++) {outFile << v[ii] << "\n";}
}

void dump_result_to_disk(const std::string base_loc,
    const std::string result_name, const double *v, const int N,
    const int index)
{
    std::ostringstream stringStream;
    stringStream << base_loc << "/" << result_name;
    const std::string base_loc_str = stringStream.str();
    const std::string path = filename(base_loc_str, index);
    dump_arr_of_doubles_to_disk(path, v, N);
}


void dump_arr_of_ints_to_disk(const std::string path,
    const int *v, const int N)
{
    std::ofstream outFile(path);
    for (int ii=0; ii<N; ii++)
    {
        outFile << v[ii] << "\n";
    }
}


void dump_result_ints_to_disk(const std::string base_loc,
    const std::string result_name, const int *v, const int N,
    const int index)
{
    std::ostringstream stringStream;
    stringStream << base_loc << "/" << result_name;
    const std::string base_loc_str = stringStream.str();
    const std::string path = filename(base_loc_str, index);
    dump_arr_of_ints_to_disk(path, v, N);
}


void dump_grid_to_disk(const std::string base_loc,
    const std::string grid_name, const double *v, const int N)
{
    std::ostringstream stringStream;
    stringStream << base_loc << "/" << grid_name << ".txt";
    const std::string path = stringStream.str();
    dump_arr_of_doubles_to_disk(path, v, N);
}
