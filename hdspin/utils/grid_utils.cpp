#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

#include "general_utils.h"
#include "grid_utils.h"


// ============================================================================
// Energy =====================================================================
// ============================================================================

EnergyGrid::EnergyGrid(const std::string grid_location)
{
    // printf("Initializing grid from %s\n", grid_location.c_str());
    std::ifstream myfile (grid_location);
    std::string line;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            grid.push_back(stoll(line));
        }
        myfile.close();
    }
    length = grid.size();
    max_time = grid[length - 1];
    // printf("Grid initialized with min/max times %lli/%lli and len %i\n",
    //     grid[0], max_time, length);
}

void EnergyGrid::open_outfile(const std::string d)
{
    outfile = fopen(d.c_str(), "w");
}

void EnergyGrid::close_outfile(){fclose(outfile);}

void EnergyGrid::step(const double current_time, const long long config_int,
    const long long config_IS_int, const double energy, const double energy_IS)
{

    // No updates necessary
    if (current_time <= grid[pointer]){return;}

    if (pointer > length - 1){return;}

    // Write to the outfile
    while (grid[pointer] < current_time)
    {   
        fprintf(outfile, "%lli %lli %.05f %lli %.05f\n", grid[pointer],
            config_int, energy, config_IS_int, energy_IS);
        pointer += 1;
        if (pointer > length - 1){return;}
    }
}



// ============================================================================
// Psi ========================================================================
// ============================================================================


PsiConfigCounter::PsiConfigCounter(const int log_N_timesteps,
    const std::string outfile_loc)
{
    max_counter = (long long) log2l(ipow(10, log_N_timesteps));

    // Give the max counter a lot of space
    max_counter += 10;

    for (int ii=0; ii<max_counter; ii++){counter.push_back(0);}
    outfile_location = outfile_loc;
}

void PsiConfigCounter::write_to_disk()
{
    FILE *outfile;
    outfile = fopen(outfile_location.c_str(), "w");
    for (int ii=0; ii<max_counter; ii++)
    {
        fprintf(outfile, "%i %lli\n", ii, counter[ii]);
    }
    fclose(outfile);
}

void PsiConfigCounter::step(const long double t)
{

    long long key;

    // If the waiting time is < 1, round it to 1.
    if (t < 1.0){key = 0;}

    else
    {
        const long double log_t = log2l(t);
        key = (long long) roundl(log_t);

        // We ignore any crazy waiting times produced near the end of the
        // Gillespie dynamics since they can be chalked up to edge effects.
        if (key > max_counter - 1){return;}
    }
    counter[key] += 1;
}
