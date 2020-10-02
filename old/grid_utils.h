#ifndef GRID_UTILS_H
#define GRID_UTILS_H

#include <iostream>
#include <vector>
#include <fstream>

#include "structure_utils.h"


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
    void step(const double, const SystemInformation, const SystemInformation);
    
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
// Psi config =================================================================
// ============================================================================


class PsiConfigCounter
{
public:

    // Constructor
    PsiConfigCounter(const int);

    // Close the output stream when done; this also writes the counter to
    // disk, since this file will be relatively small compared to the others
    // that are written dynamically.
    void write_to_disk(const std::string);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step(const long double, const bool);
    
private:

    // The log2-binned trapping times
    std::vector<long long> counter, counter_IS;

    // The maximum key value
    long long max_counter;

    // The location to save the file
    std::string outfile_location;
};



// ============================================================================
// Pi/Aging config ============================================================
// ============================================================================

class AgingConfigGrid
{
public:

    // Constructor: reads in both the pi1 and pi2 grids from disk.
    AgingConfigGrid(const std::string);

    // Want to open the filestream to the save directory
    void open_outfile(const std::string, const std::string);

    // Close the output stream when done
    void close_outfile();

    void step(const double, const long long, const long long,
        const long long);

private:

    // The location to save the file
    std::string outfile_location;

    // The pi 1 and 2 grids
    std::vector<long long> grid_pi1;
    std::vector<long long> grid_pi2;
    int length_pi1, length_pi2;
    long long max_time_pi1, max_time_pi2;

    // The outstream for this tracker
    FILE *outfile_pi1;
    FILE *outfile_pi2;

    // Define the pointers
    int pointer1 = 0;
    int pointer2 = 0;
};


// ============================================================================
// Psi basin =================================================================
// ============================================================================


class PsiBasinCounter
{
public:

    // Constructor
    PsiBasinCounter(const int, const double, const double);

    // Close the output stream when done; this also writes the counter to
    // disk, since this file will be relatively small compared to the others
    // that are written dynamically.
    void write_to_disk(const std::string);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step(SystemInformation *, SystemInformation *);
    
private:

    // The thresholds
    double energetic_threshold, entropic_attractor;

    // The log2-binned trapping times
    std::vector<long long> counter_E, counter_E_IS, counter_S, counter_S_IS;

    // The maximum key value
    long long max_counter;

    // The location to save the file
    std::string outfile_location;
};



#endif
