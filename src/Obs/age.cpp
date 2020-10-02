#include "Obs/base.h"
#include "Obs/age.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

AgingConfig::AgingConfig(const FileNames fnames)
{
    const std::string pi_1_grid_location = fnames.grids_directory + "/pi1.txt";
    const std::string pi_2_grid_location = fnames.grids_directory + "/pi2.txt";

    load_long_long_grid_(grid_pi1, pi_1_grid_location);
    length_pi1 = grid_pi1.size();
    max_time_pi1 = grid_pi1[length_pi1 - 1];

    load_long_long_grid_(grid_pi2, pi_2_grid_location);
    length_pi2 = grid_pi2.size();
    max_time_pi2 = grid_pi2[length_pi2 - 1];

    outfile_pi1 = fopen(fnames.aging_config_1.c_str(), "w");
    outfile_pi2 = fopen(fnames.aging_config_2.c_str(), "w");
}


void AgingConfig::_help_step_1_(const long double simulation_clock,
    const long long config_index, const Vals prev)
{
    while (grid_pi1[pointer1] < simulation_clock)
    {   
        fprintf(outfile_pi1, "%lli %lli %lli %lli\n",
            grid_pi1[pointer1], config_index, prev.int_rep, prev.int_rep_IS);
        pointer1 += 1;
        if (pointer1 > length_pi1 - 1){break;}
    }  
}


void AgingConfig::_help_step_2_(const long double simulation_clock,
    const long long config_index, const Vals prev)
{
    while (grid_pi2[pointer2] < simulation_clock)
    {   
        fprintf(outfile_pi2, "%lli %lli %lli %lli\n",
            grid_pi2[pointer2], config_index, prev.int_rep, prev.int_rep_IS);
        pointer2 += 1;
        if (pointer2 > length_pi2 - 1){break;}
    }  
}


// Note that the config index is iterated in the main simulation loop every
// time a new configuration is accepted.
void AgingConfig::step_(const long double simulation_clock,
    const long long config_index, const Vals prev)
{
    if (simulation_clock <= grid_pi1[pointer1]){;}
    else if (pointer1 > length_pi1 - 1){;}
    else{_help_step_1_(simulation_clock, config_index, prev);}

    if (simulation_clock <= grid_pi2[pointer2]){;}
    else if (pointer2 > length_pi2 - 1){;}
    else{_help_step_2_(simulation_clock, config_index, prev);}
}

AgingConfig::~AgingConfig()
{
    fclose(outfile_pi1);
    fclose(outfile_pi2);
}
