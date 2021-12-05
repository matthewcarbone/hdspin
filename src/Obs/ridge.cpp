#include "Obs/ridge.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

RidgeBase::RidgeBase(const FileNames fnames,
    const RuntimeParameters rtp) : rtp(rtp) {}


void RidgeBase::_log_ridge(const double current_energy)
{

    if (_logged >= rtp.max_ridges){return;}

    const int diff_status = int(current_energy == _last_energy);

    fprintf(outfile, "%.05e %i %0.5e %i\n",
        _current_ridge, _steps_above, _time_above, diff_status);

    _logged += 1;
}


void RidgeBase::step(const Vals prev, const Vals curr, const double waiting_time)
{

    if (!_threshold_valid){return;}

    const double _prev_energy = prev.energy_IS;
    const double _curr_energy = curr.energy_IS;
    const double _prev_energy_compare = prev.energy_IS;
    const double _curr_energy_compare = curr.energy_IS;

    // Just exited a basin, track the previous energy
    if ((_prev_energy < _threshold) && (_curr_energy >= _threshold))
    {
        // Initialize the last energy before the basin was exited via one of
        // two ways. If inherent_structure == 0, we set the last energy to that
        // of the previous in the standard trajectory. If inherent_structure
        // != 0, then the last energy should be that of the inherent structure.
        // This is true even when inherent_structure == 2, which takes the
        // _prev_energy (for comparing to the threshold) to be the standard
        // trajectory, but sets the last energy according to the IS.
        _last_energy = _prev_energy_compare;
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
        if (_exited_first_basin){_log_ridge(_curr_energy_compare);}
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

