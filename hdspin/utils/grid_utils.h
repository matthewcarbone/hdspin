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

    // Get's the grid
    std::vector<long long> get_grid();

    // Want to open the filestream to the save directory
    void open_outfile(const std::string);

    // Close the output stream when done
    void close_outfile();

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step(const double, const double, const int *, const int,
        const double *, long long *);
    
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

#endif

