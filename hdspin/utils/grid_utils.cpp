#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

#include "general_utils.h"
#include "grid_utils.h"
#include "structure_utils.h"


// ============================================================================
// Energy =====================================================================
// ============================================================================

EnergyGrid::EnergyGrid(const std::string grid_directory)
{
    // printf("Initializing grid from %s\n", grid_location.c_str());
    const std::string grid_location = grid_directory + "/energy.txt";
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

void EnergyGrid::step(const double current_time, const SystemInformation sys,
    const SystemInformation inh)
{

    // No updates necessary
    if (current_time <= grid[pointer]){return;}

    if (pointer > length - 1){return;}

    // Write to the outfile
    while (grid[pointer] < current_time)
    {   
        fprintf(outfile, "%lli %lli %.05Lf %lli %.05Lf\n", grid[pointer],
            sys.x, sys.e, inh.x, inh.e);
        pointer += 1;
        if (pointer > length - 1){return;}
    }
}



// ============================================================================
// Psi config =================================================================
// ============================================================================


PsiConfigCounter::PsiConfigCounter(const int log_N_timesteps)
{
    max_counter = (long long) log2l(ipow(10, log_N_timesteps));

    // Give the max counter a lot of space
    max_counter += 10;

    for (int ii=0; ii<max_counter; ii++){counter.push_back(0);}
    for (int ii=0; ii<max_counter; ii++){counter_IS.push_back(0);}
}

void PsiConfigCounter::write_to_disk(const std::string outfile_location)
{
    FILE *outfile;
    outfile = fopen(outfile_location.c_str(), "w");
    for (int ii=0; ii<max_counter; ii++)
    {
        fprintf(outfile, "%i %lli %lli\n", ii, counter[ii], counter_IS[ii]);
    }
    fclose(outfile);
}

void PsiConfigCounter::step(const long double t, const bool inherent_structure)
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

    if (inherent_structure == true){counter_IS[key] += 1;}
    else{counter[key] += 1;}
}


// ============================================================================
// Pi/Aging config ============================================================
// ============================================================================


AgingConfigGrid::AgingConfigGrid(const std::string grid_directory)
{
    const std::string pi_1_grid_location = grid_directory + "/pi1.txt";
    const std::string pi_2_grid_location = grid_directory + "/pi2.txt";

    std::string line1;

    std::ifstream myfile_1 (pi_1_grid_location);
    if (myfile_1.is_open())
    {
        while (getline(myfile_1, line1))
        {
            grid_pi1.push_back(stoll(line1));
        }
        myfile_1.close();
    }
    length_pi1 = grid_pi1.size();
    max_time_pi1 = grid_pi1[length_pi1 - 1];

    std::string line2;

    std::ifstream myfile_2 (pi_2_grid_location);
    if (myfile_2.is_open())
    {
        while (getline(myfile_2, line2))
        {
            grid_pi2.push_back(stoll(line2));
        }
        myfile_2.close();
    }
    length_pi2 = grid_pi2.size();
    max_time_pi2 = grid_pi2[length_pi2 - 1];
}

void AgingConfigGrid::open_outfile(const std::string pi1_path,
    const std::string pi2_path)
{
    outfile_pi1 = fopen(pi1_path.c_str(), "w");
    outfile_pi2 = fopen(pi2_path.c_str(), "w");
}

void AgingConfigGrid::close_outfile()
{
    fclose(outfile_pi1);
    fclose(outfile_pi2);
}


void AgingConfigGrid::step(const double current_time,
    const long long config_index, const long long config_int,
    const long long config_IS_int)
{
    if (current_time <= grid_pi1[pointer1]){;}
    else if (pointer1 > length_pi1 - 1){;}
    else
    {
        // Write to the outfile
        while (grid_pi1[pointer1] < current_time)
        {   
            fprintf(outfile_pi1, "%lli %lli %lli %lli\n",
                grid_pi1[pointer1], config_index, config_int, config_IS_int);
            pointer1 += 1;
            if (pointer1 > length_pi1 - 1){break;}
        }  
    }

    if (current_time <= grid_pi2[pointer2]){;}
    else if (pointer2 > length_pi2 - 1){;}
    else
    {
        // Write to the outfile
        while (grid_pi2[pointer2] < current_time)
        {   
            fprintf(outfile_pi2, "%lli %lli %lli %lli\n",
                grid_pi2[pointer2], config_index, config_int, config_IS_int);
            pointer2 += 1;
            if (pointer2 > length_pi2 - 1){break;}
        }  
    }
}
