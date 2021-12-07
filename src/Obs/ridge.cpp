#include "Obs/ridge.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

RidgeBase::RidgeBase(const FileNames fnames,
    const RuntimeParameters rtp) : rtp(rtp), fnames(fnames) {}


void RidgeBase::_log_ridge(const double current_energy)
{

    const int diff_status = int(current_energy == _last_energy);

    fprintf(outfile, "%.05e %i %0.5e %i\n",
        _current_ridge, _steps_above, _time_above, diff_status);

    _logged += 1;
}


void RidgeBase::step(const Vals prev, const Vals curr, const double waiting_time)
{

    if (!_threshold_valid){return;}
    if (_logged >= rtp.max_ridges){return;}

    const double _prev_energy = prev.energy;
    const double _curr_energy = curr.energy;

    // Just exited a basin, track the previous energy
    if ((_prev_energy < _threshold) && (_curr_energy >= _threshold))
    {
        // The energy before exiting the basin is defined as the inherent
        // structure energy if we are using a system with memory. Else, we just
        // use the energy itself
        _last_energy = (rtp.memory != 0) ? prev.energy_IS : prev.energy;
        _current_ridge = _curr_energy;
        _exited_first_basin = true;
    }

    // Still above the basin. The current ridge energy is defined as the
    // maximum energy reached above the threshold.
    else if ((_prev_energy >= _threshold) && (_curr_energy >= _threshold))
    {
        _current_ridge =
            _curr_energy > _current_ridge ? _curr_energy : _current_ridge;
        _steps_above += 1;
        _time_above += waiting_time;
    }

    // Just dropped back below the threshold: log the current ridge energy
    else if ((_prev_energy >= _threshold) && (_curr_energy < _threshold))
    {
        _steps_above += 1;
        _time_above += waiting_time;
        if (_exited_first_basin)
        {
            const double _curr_energy_compare =
                (rtp.memory != 0) ? curr.energy_IS : curr.energy;
            _log_ridge(_curr_energy_compare);
        }
        _steps_above = 0;
        _time_above = 0.0;
    }

}

RidgeBase::~RidgeBase()
{
    fclose(outfile);
}


RidgeEnergy::RidgeEnergy(const FileNames fnames,
    const RuntimeParameters rtp) : RidgeBase(fnames, rtp)
{
    outfile = fopen(fnames.ridge_E_all.c_str(), "w");
    _threshold = rtp.energetic_threshold;
}


RidgeAttractor::RidgeAttractor(const FileNames fnames,
    const RuntimeParameters rtp) : RidgeBase(fnames, rtp)
{
    outfile = fopen(fnames.ridge_S_all.c_str(), "w");
    _threshold = rtp.entropic_attractor;
    _threshold_valid = rtp.valid_entropic_attractor;
}

