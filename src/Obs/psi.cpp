#include <set>

#include "Obs/psi.h"
#include "Utils/utils.h"
#include "Utils/structures.h"


long long _get_key(const long double local_waiting_time)
{
    long long key;

    // If the waiting time is <= 1, round it to 1.
    if (local_waiting_time <= 1.0){key = 0;}

    else
    {
        const long double log_t = log2l(local_waiting_time);
        key = (long long) roundl(log_t);
    }

    return key;
}


PsiConfigBase::PsiConfigBase(const FileNames fnames,
    const RuntimeParameters rtp) : fnames(fnames), rtp(rtp)
{
    _max_counter = (long long) log2l(ipow(10, rtp.log_N_timesteps));

    // Give the max counter a lot of space
    _max_counter += 10;

    for (int ii=0; ii<_max_counter; ii++){_counter.push_back(0);}
}


void PsiConfigBase::_help_step()
{
    const long long key = _get_key(_waiting_time);
    _waiting_time = 0.0;

    // We ignore any crazy waiting times produced near the end of the
    // Gillespie dynamics since they can be chalked up to edge effects.
    if (key > _max_counter - 1){return;}

    // Update
    _counter[key] += 1;
}


// Psi Config -----------------------------------------------------------------

PsiConfig::PsiConfig(const FileNames fnames,
    const RuntimeParameters rtp) : PsiConfigBase(fnames, rtp){}

void PsiConfig::step(const long double current_waiting_time, const Vals prev,
    const Vals curr)
{
    // No matter what, we update the internal states. If the configs are logged
    // to the counters they are reset in the helper functions.
    _waiting_time += current_waiting_time;
    if (curr.int_rep != prev.int_rep){_help_step();}
}

PsiConfig::~PsiConfig()
{
    outfile = fopen(fnames.psi_config.c_str(), "w");
    for (int ii=0; ii<_max_counter; ii++)
    {
        fprintf(outfile, "%lli\n", _counter[ii]);
    }
    fclose(outfile);
}


// Psi Config Inherent Structure ----------------------------------------------

PsiConfigInherentStructure::PsiConfigInherentStructure(const FileNames fnames,
    const RuntimeParameters rtp) : PsiConfigBase(fnames, rtp){}

void PsiConfigInherentStructure::step(
    const long double current_waiting_time, const Vals prev, const Vals curr)
{
    if (rtp.memory == 0){return;}

    // No matter what, we update the internal states. If the configs are logged
    // to the counters they are reset in the helper functions.
    _waiting_time += current_waiting_time;
    if (curr.int_rep_IS != prev.int_rep_IS){_help_step();}
}

PsiConfigInherentStructure::~PsiConfigInherentStructure()
{
    // Don't do anything if we're using memoryless dynamics
    if (rtp.memory == 0){return;}
    outfile = fopen(fnames.psi_config_IS.c_str(), "w");
    for (int ii=0; ii<_max_counter; ii++)
    {
        fprintf(outfile, "%lli\n", _counter[ii]);
    }
    fclose(outfile);
}


// Psi Basin ------------------------------------------------------------------

PsiBasinBase::PsiBasinBase(const FileNames fnames,
    const RuntimeParameters rtp) : fnames(fnames), rtp(rtp)
{
    _max_counter = (long long) log2l(ipow(10, rtp.log_N_timesteps));

    // Give the max counter a lot of space
    _max_counter += 10;

    for (int ii=0; ii<_max_counter; ii++){_counter.push_back(0);}
    for (int ii=0; ii<_max_counter; ii++)
    {
        _counter_unique_configs_per_basin.push_back(0);
    }
}


void PsiBasinBase::_help_step(const double _prev_energy,
    const double _curr_energy, const long long _prev_int_rep,
    const long double current_waiting_time)
{
    if (_prev_energy < _threshold)
    {
        _waiting_time += current_waiting_time;
        _tmp_unique_configs_in_basin.push_back(_prev_int_rep);

        if (_curr_energy >= _threshold)
        {
            const long double local_waiting_time = _waiting_time;
            _waiting_time = 0.0;

            const long long key = _get_key(local_waiting_time);
            if (key <= _max_counter - 1){_counter[key] += 1;}

            // We also count the number of unique configs per basin
            const long double n_unique =
                std::set<long long>(_tmp_unique_configs_in_basin.begin(),
                    _tmp_unique_configs_in_basin.end()).size();
            const long long key2 = _get_key(n_unique);
            if (key2 <= _max_counter - 1)
            {
                _counter_unique_configs_per_basin[key2] += 1;
            }
            _tmp_unique_configs_in_basin.clear();
        }
    }
}

// Default step function
void PsiBasinBase::step(const long double current_waiting_time,
    const Vals prev, const Vals curr)
{
    if (!_threshold_valid){return;}
    _help_step(prev.energy, curr.energy, prev.int_rep, current_waiting_time);
}


// Dump the outfile
void PsiBasinBase::_dump_outfile()
{
    for (int ii=0; ii<_max_counter; ii++)
    {
        fprintf(outfile, "%lli %lli\n",
            _counter[ii], _counter_unique_configs_per_basin[ii]);
    }
    fclose(outfile);
}


// Psi Basin (E) --------------------------------------------------------------

PsiBasinThreshold::PsiBasinThreshold(const FileNames fnames,
    const RuntimeParameters rtp) : PsiBasinBase(fnames, rtp)
{
    _threshold = rtp.energetic_threshold;
}


PsiBasinThreshold::~PsiBasinThreshold()
{
    outfile = fopen(fnames.psi_basin_E.c_str(), "w");
    _dump_outfile();
}


// Psi Basin (S) --------------------------------------------------------------

PsiBasinAttractor::PsiBasinAttractor(const FileNames fnames,
    const RuntimeParameters rtp) : PsiBasinBase(fnames, rtp)
{
    _threshold = rtp.entropic_attractor;
}


PsiBasinAttractor::~PsiBasinAttractor()
{
    if (!_threshold_valid){return;}
    outfile = fopen(fnames.psi_basin_S.c_str(), "w");
    _dump_outfile();
}


// Psi Basin (IS+E) -----------------------------------------------------------

PsiBasinThresholdInherentStructure::PsiBasinThresholdInherentStructure(
    const FileNames fnames, const RuntimeParameters rtp) : PsiBasinBase(fnames, rtp)
{
    _threshold = rtp.energetic_threshold;
}

void PsiBasinThresholdInherentStructure::step(
    const long double current_waiting_time, const Vals prev, const Vals curr)
{
    if (rtp.memory == 0){return;}
    _help_step(prev.energy_IS, curr.energy_IS, prev.int_rep_IS,
        current_waiting_time);
}

PsiBasinThresholdInherentStructure::~PsiBasinThresholdInherentStructure()
{
    if (rtp.memory == 0){return;}
    outfile = fopen(fnames.psi_basin_E_IS.c_str(), "w");
    _dump_outfile();
}


// Psi Basin (IS+S) -----------------------------------------------------------

PsiBasinAttractorInherentStructure::PsiBasinAttractorInherentStructure(
    const FileNames fnames, const RuntimeParameters rtp) : PsiBasinBase(fnames, rtp)
{
    _threshold = rtp.entropic_attractor;
}

void PsiBasinAttractorInherentStructure::step(
    const long double current_waiting_time, const Vals prev, const Vals curr)
{
    if (!_threshold_valid){return;}
    if (rtp.memory == 0){return;}
    _help_step(prev.energy_IS, curr.energy_IS, prev.int_rep_IS,
        current_waiting_time);
}

PsiBasinAttractorInherentStructure::~PsiBasinAttractorInherentStructure()
{
    if (!_threshold_valid){return;}
    if (rtp.memory == 0){return;}
    outfile = fopen(fnames.psi_basin_E_IS.c_str(), "w");
    _dump_outfile();
}



