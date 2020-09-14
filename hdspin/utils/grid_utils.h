#include <iostream>
#include <vector>
#include <fstream>

#ifndef GRID_UTILS_H
#define GRID_UTILS_H


// ============================================================================
// Energy =====================================================================
// ============================================================================

class EnergyGrid
{
public:

    // Constructor: reads in the grid
    EnergyGrid(const std::string);

    // Want to open the filestream to the save directory
    void open_outfile(const std::string);

    // Close the output stream when done
    void close_outfile();

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step(const double, const long long, const long long, const double,
        const double);
    
private:

    // Points to the last-saved point on the grid
    int pointer = 0;

    // The maximum time on the grid
    long long max_time;

    // The grid (in time) itself
    std::vector<long long> grid;

    // The outstream for this tracker
    FILE *outfile;

    // The length of the grid
    int length;
};




// ============================================================================
// Psi ========================================================================
// ============================================================================


class PsiConfigCounter
{
public:

    // Constructor: reads in the grid
    PsiConfigCounter(const int, const std::string outfile_loc);

    // Close the output stream when done; this also writes the counter to
    // disk, since this file will be relatively small compared to the others
    // that are written dynamically.
    void write_to_disk();

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step(const long double);
    
private:

    // The log2-binned trapping times
    std::vector<long long> counter;

    // The maximum key value
    long long max_counter;

    // The location to save the file
    std::string outfile_location;
};




#endif
