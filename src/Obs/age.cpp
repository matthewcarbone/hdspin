#include <cassert>

#include "Obs/age.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

Aging::Aging(const FileNames fnames,
    const RuntimeParameters rtp) : fnames(fnames), rtp(rtp)
{
    const std::string pi_1_grid_location = fnames.grids_directory + "/pi1.txt";
    const std::string pi_2_grid_location = fnames.grids_directory + "/pi2.txt";

    load_long_long_grid_(grid_pi1, pi_1_grid_location);
    length = grid_pi1.size();

    load_long_long_grid_(grid_pi2, pi_2_grid_location);
    const int length_pi2 = grid_pi2.size();

    assert(length == length_pi2);
}

void Aging::_help_step_1(const long double simulation_clock,
    const long long int_rep)
{
    while (grid_pi1[pointer1] < simulation_clock)
    {
        results1.push_back(int_rep);
        pointer1 += 1;
        if (pointer1 > length - 1){break;}
    }
}


void Aging::_help_step_2(const long double simulation_clock,
    const long long int_rep)
{
    while (grid_pi2[pointer2] < simulation_clock)
    {
        results2.push_back(int_rep);
        pointer2 += 1;
        if (pointer2 > length - 1){break;}
    }
}


void Aging::_dump_outfile()
{
    for (int ii = 0; ii < length; ii++)
    {
        fprintf(outfile, "%lli %lli\n", results1[ii], results2[ii]);
    }

    fclose(outfile);
}



// Aging config ---------------------------------------------------------------

AgingConfig::AgingConfig(const FileNames fnames,
    const RuntimeParameters rtp) : Aging(fnames, rtp)
{
    outfile = fopen(fnames.aging_config.c_str(), "w");
}


void AgingConfig::step(const long double simulation_clock, const Vals prev)
{
    if (simulation_clock <= grid_pi1[pointer1]){;}
    else if (pointer1 > length - 1){;}
    else{_help_step_1(simulation_clock, prev.int_rep);}

    if (simulation_clock <= grid_pi2[pointer2]){;}
    else if (pointer2 > length - 1){;}
    else{_help_step_2(simulation_clock, prev.int_rep);}
}


AgingConfig::~AgingConfig()
{
    _dump_outfile();
}



// Aging config IS ------------------------------------------------------------

AgingConfigInherentStructure::AgingConfigInherentStructure(const FileNames fnames,
    const RuntimeParameters rtp) : Aging(fnames, rtp)
{
    if (rtp.memory != 0){outfile = fopen(fnames.aging_config_IS.c_str(), "w");}
}


void AgingConfigInherentStructure::step(const long double simulation_clock,
    const Vals prev)
{
    if (rtp.memory == 0){return;}

    if (simulation_clock <= grid_pi1[pointer1]){;}
    else if (pointer1 > length - 1){;}
    else{_help_step_1(simulation_clock, prev.int_rep_IS);}

    if (simulation_clock <= grid_pi2[pointer2]){;}
    else if (pointer2 > length - 1){;}
    else{_help_step_2(simulation_clock, prev.int_rep_IS);}
}


AgingConfigInherentStructure::~AgingConfigInherentStructure()
{
    if (rtp.memory != 0){_dump_outfile();}
}




void AgingBasin::_dump_outfile()
{
    for (int ii = 0; ii < length; ii++)
    {
        fprintf(outfile, "%lli %i %lli %i\n", vec_basin_index_1[ii],
            vec_prev_state_in_basin_1[ii], vec_basin_index_2[ii],
            vec_prev_state_in_basin_2[ii]);
    }
    fclose(outfile);
}

AgingBasin::AgingBasin(const FileNames fnames, const RuntimeParameters rtp) :
    Aging(fnames, rtp) {}


void AgingBasin::_help_step_1_(const long double simulation_clock,
    const double prev_E)
{
    const int prev_state_in_basin = prev_E < _threshold ? 1 : 0;

    // Column 1: basin index
    // Column 2: in basin or not (1 if in basin, 0 if not in basin)
    while (grid_pi1[pointer1] < simulation_clock)
    {
        vec_basin_index_1.push_back(basin_index_1);
        vec_prev_state_in_basin_1.push_back(prev_state_in_basin);
        pointer1 += 1;
        if (pointer1 > length - 1){break;}
    }
}


void AgingBasin::_help_step_2_(const long double simulation_clock,
    const double prev_E)
{
    const int prev_state_in_basin = prev_E < _threshold ? 1 : 0;

    while (grid_pi2[pointer2] < simulation_clock)
    {
        vec_basin_index_2.push_back(basin_index_2);
        vec_prev_state_in_basin_2.push_back(prev_state_in_basin);
        pointer2 += 1;
        if (pointer2 > length - 1){break;}
    }
}


void AgingBasin::_help_step(const long double simulation_clock,
    const double prev_E, const double curr_E)
{
    if (pointer1 > length - 1){;}
    else
    {
        if (simulation_clock > grid_pi1[pointer1])
        {
            _help_step_1_(simulation_clock, prev_E);
        }
        // Just entered a basin, iterate the basin index. We do this after
        // stepping the grids because the tracer is assumed to be in this
        // configuration until, but not including, the current time indexed by
        //the simulation clock.
        if ((prev_E >= _threshold) && (curr_E < _threshold))
        {
            basin_index_1++;
        }
    }

    if (pointer2 > length - 1){;}
    else
    {
        if (simulation_clock > grid_pi2[pointer2])
        {
            _help_step_2_(simulation_clock, prev_E);
        }

        if ((prev_E >= _threshold) && (curr_E < _threshold))
        {
            basin_index_2++;
        }
    }
}


void AgingBasin::step(const long double simulation_clock,
    const Vals prev, const Vals curr)
{
    _help_step(simulation_clock, prev.energy, curr.energy);
}



AgingBasinThreshold::AgingBasinThreshold(const FileNames fnames,
    const RuntimeParameters rtp) : AgingBasin(fnames, rtp)
{
    outfile = fopen(fnames.aging_basin_E.c_str(), "w");
    _threshold = rtp.energetic_threshold;
};

AgingBasinThreshold::~AgingBasinThreshold()
{
    _dump_outfile();
}


///////


AgingBasinAttractor::AgingBasinAttractor(const FileNames fnames,
    const RuntimeParameters rtp) : AgingBasin(fnames, rtp)
{
    if (!rtp.valid_entropic_attractor){return;}
    outfile = fopen(fnames.aging_basin_S.c_str(), "w");
    _threshold = rtp.entropic_attractor;
};

void AgingBasinAttractor::step(const long double simulation_clock,
    const Vals prev, const Vals curr)
{
    if (!rtp.valid_entropic_attractor){return;}
    _help_step(simulation_clock, prev.energy, curr.energy);
}

AgingBasinAttractor::~AgingBasinAttractor()
{
    if (!rtp.valid_entropic_attractor){return;}
    _dump_outfile();
}


///////


AgingBasinThresholdInherentStructure::AgingBasinThresholdInherentStructure(
    const FileNames fnames, const RuntimeParameters rtp) : AgingBasin(fnames, rtp)
{
    if (rtp.memory == 0){return;}
    outfile = fopen(fnames.aging_basin_E_IS.c_str(), "w");
    _threshold = rtp.energetic_threshold;
};

void AgingBasinThresholdInherentStructure::step(
    const long double simulation_clock, const Vals prev, const Vals curr)
{
    if (rtp.memory == 0){return;}
    _help_step(simulation_clock, prev.energy_IS, curr.energy_IS);
}

AgingBasinThresholdInherentStructure::~AgingBasinThresholdInherentStructure()
{
    if (rtp.memory == 0){return;}
    _dump_outfile();
}


///////


AgingBasinAttractorInherentStructure::AgingBasinAttractorInherentStructure(
    const FileNames fnames, const RuntimeParameters rtp) : AgingBasin(fnames, rtp)
{
    if (rtp.memory == 0 || !rtp.valid_entropic_attractor){return;}
    outfile = fopen(fnames.aging_basin_S_IS.c_str(), "w");
    _threshold = rtp.entropic_attractor;
};

void AgingBasinAttractorInherentStructure::step(const long double simulation_clock,
    const Vals prev, const Vals curr)
{
    if (rtp.memory == 0 || !rtp.valid_entropic_attractor){return;}
    _help_step(simulation_clock, prev.energy_IS, curr.energy_IS);
}

AgingBasinAttractorInherentStructure::~AgingBasinAttractorInherentStructure()
{
    if (rtp.memory == 0 || !rtp.valid_entropic_attractor){return;}
    _dump_outfile();
}

