#include <cassert>

#include "Obs/base.h"
#include "Obs/age.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

Aging::Aging(const FileNames fnames)
{
    const std::string pi_1_grid_location = fnames.grids_directory + "/pi1.txt";
    const std::string pi_2_grid_location = fnames.grids_directory + "/pi2.txt";

    load_long_long_grid_(grid_pi1, pi_1_grid_location);
    length_pi1 = grid_pi1.size();
    max_time_pi1 = grid_pi1[length_pi1 - 1];

    load_long_long_grid_(grid_pi2, pi_2_grid_location);
    length_pi2 = grid_pi2.size();
    max_time_pi2 = grid_pi2[length_pi2 - 1];

    assert(length_pi1 == length_pi2);
}


Aging::~Aging()
{
    fclose(outfile_pi1);
    fclose(outfile_pi2);
}


AgingConfig::AgingConfig(const FileNames fnames) : Aging(fnames)
{
    outfile_pi1 = fopen(fnames.aging_config_1.c_str(), "w");
    outfile_pi2 = fopen(fnames.aging_config_2.c_str(), "w");
};


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


AgingBasin::AgingBasin(const FileNames fnames, const RuntimeParameters rtp) :
    Aging(fnames), rtp(rtp)
{
    outfile_pi1 = fopen(fnames.aging_basin_1.c_str(), "w");
    outfile_pi2 = fopen(fnames.aging_basin_2.c_str(), "w");  
};



void AgingBasin::_help_step_1_(const long double simulation_clock,
    const Vals prev, const Vals curr)
{
    const int prev_state_in_basin_E =
        prev.energy < rtp.energetic_threshold ? 1 : 0;
    const int prev_IS_state_in_basin_E =
        prev.energy_IS < rtp.energetic_threshold ? 1 : 0;
    const int prev_state_in_basin_S =
        prev.energy < rtp.entropic_attractor ? 1 : 0;
    const int prev_IS_state_in_basin_S =
        prev.energy_IS < rtp.entropic_attractor ? 1 : 0;

    // Column 1: grid point
    // ENERGY:
    // Column 2: basin index
    // Column 3: in basin or not (1 if in basin, 0 if not in basin)
    // Column 4: IS basin index
    // Column 5: IS in basin or not
    // ENTROPY:
    // Column 6: basin index
    // Column 7: in basin or not (1 if in basin, 0 if not in basin)
    // Column 8: IS basin index
    // Column 9: IS in basin or not
    while (grid_pi1[pointer1] < simulation_clock)
    {   
        fprintf(outfile_pi1, "%lli %lli %i %lli %i %lli %i %lli %i\n",
            grid_pi1[pointer1], bi1_E, prev_state_in_basin_E, bi1_E_IS,
            prev_IS_state_in_basin_E, bi1_S, prev_state_in_basin_S, bi1_S_IS,
            prev_IS_state_in_basin_S);
        pointer1 += 1;
        if (pointer1 > length_pi1 - 1){break;}
    }
}


void AgingBasin::_help_step_2_(const long double simulation_clock,
    const Vals prev, const Vals curr)
{
    const int prev_state_in_basin_E =
        prev.energy < rtp.energetic_threshold ? 1 : 0;
    const int prev_IS_state_in_basin_E =
        prev.energy_IS < rtp.energetic_threshold ? 1 : 0;
    const int prev_state_in_basin_S =
        prev.energy < rtp.entropic_attractor ? 1 : 0;
    const int prev_IS_state_in_basin_S =
        prev.energy_IS < rtp.entropic_attractor ? 1 : 0;

    // Column 1: grid point
    // ENERGY:
    // Column 2: basin index
    // Column 3: in basin or not (1 if in basin, 0 if not in basin)
    // Column 4: IS basin index
    // Column 5: IS in basin or not
    // ENTROPY:
    // Column 6: basin index
    // Column 7: in basin or not (1 if in basin, 0 if not in basin)
    // Column 8: IS basin index
    // Column 9: IS in basin or not
    while (grid_pi2[pointer2] < simulation_clock)
    {   
        fprintf(outfile_pi2, "%lli %lli %i %lli %i %lli %i %lli %i\n",
            grid_pi2[pointer2], bi2_E, prev_state_in_basin_E, bi2_E_IS,
            prev_IS_state_in_basin_E, bi2_S, prev_state_in_basin_S, bi2_S_IS,
            prev_IS_state_in_basin_S);
        pointer2 += 1;
        if (pointer2 > length_pi2 - 1){break;}
    }   
}


void AgingBasin::step_(const long double simulation_clock,
    const Vals prev, const Vals curr)
{
    if (pointer1 > length_pi1 - 1){;}
    else
    {
        if (simulation_clock > grid_pi1[pointer1])
        {
            _help_step_1_(simulation_clock, prev, curr);
        }
        // Just entered a basin, iterate the basin index. We do this after
        // stepping the grids because the tracer is assumed to be in this
        // configuration until, but not including, the current time indexed by
        //the simulation clock.
        if ((prev.energy >= rtp.energetic_threshold) &&
            (curr.energy < rtp.energetic_threshold)){bi1_E++;}
        if ((prev.energy_IS >= rtp.energetic_threshold) &&
            (curr.energy_IS < rtp.energetic_threshold)){bi1_E_IS++;}
        if ((prev.energy >= rtp.entropic_attractor) &&
            (curr.energy < rtp.entropic_attractor)){bi1_S++;}
        if ((prev.energy_IS >= rtp.entropic_attractor) &&
            (curr.energy_IS < rtp.entropic_attractor)){bi1_S_IS++;}
    }

    if (pointer2 > length_pi2 - 1){;}
    else
    {
        if (simulation_clock > grid_pi2[pointer2])
        {
            _help_step_2_(simulation_clock, prev, curr);
        }
        if ((prev.energy >= rtp.energetic_threshold) &&
            (curr.energy < rtp.energetic_threshold)){bi2_E++;}
        if ((prev.energy_IS >= rtp.energetic_threshold) &&
            (curr.energy_IS < rtp.energetic_threshold)){bi2_E_IS++;}
        if ((prev.energy >= rtp.entropic_attractor) &&
            (curr.energy < rtp.entropic_attractor)){bi2_S++;}
        if ((prev.energy_IS >= rtp.entropic_attractor) &&
            (curr.energy_IS < rtp.entropic_attractor)){bi2_S_IS++;}
    }
}
