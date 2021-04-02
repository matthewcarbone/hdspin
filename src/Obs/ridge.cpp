#include <algorithm>    // std::min

#include "Obs/base.h"
#include "Obs/ridge.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

RidgeEnergy::RidgeEnergy(const FileNames fnames, const RuntimeParameters rtp,
    const bool inherent_structure, const bool energetic_threshold) :
    Base(fnames), inherent_structure(inherent_structure),
    energetic_threshold(energetic_threshold), rtp(rtp)
{
    if (energetic_threshold){threshold = rtp.energetic_threshold;}
    else{threshold = rtp.entropic_attractor;}
}


void RidgeEnergy::_log_ridge_(const double current_energy)
{
    // Logs the ridge energy

    if (current_energy == last_energy)
    {
        same.mu1 = iterative_mean(same.mu0, current_ridge, same.counter);
        same.S1 = iterative_S(same.mu0, same.mu1, current_ridge, same.S0);
        same.mu0 = same.mu1;
        same.S0 = same.S1;
        same.current_min = std::min(same.current_min, current_ridge);
        same.current_max = std::max(same.current_max, current_ridge);
        same.counter += 1;
    }
    else
    {
        diff.mu1 = iterative_mean(diff.mu0, current_ridge, diff.counter);
        diff.S1 = iterative_S(diff.mu0, diff.mu1, current_ridge, diff.S0);
        diff.mu0 = diff.mu1;
        diff.S0 = diff.S1;
        diff.current_min = std::min(diff.current_min, current_ridge);
        diff.current_max = std::max(diff.current_max, current_ridge);
        diff.counter += 1;
    }
}


void RidgeEnergy::step_(const Vals prev, const Vals curr)
{

    double _prev_energy, _curr_energy;
    if (inherent_structure)
    {
        _prev_energy = prev.energy_IS;
        _curr_energy = curr.energy_IS;
    }
    else
    {
        _prev_energy = prev.energy;
        _curr_energy = curr.energy;
    }

    // Just exited a basin, track the previous energy
    if ((_prev_energy < threshold) && (_curr_energy >= threshold))
    {
        last_energy = _prev_energy;
        current_ridge = _curr_energy;
        exited_first_basin = true;
    }

    // Still above the basin. The current ridge energy is defined as the
    // maximum energy reached above the threshold.
    else if ((_prev_energy >= threshold) && (_curr_energy >= threshold))
    {
        current_ridge = std::max(_curr_energy, current_ridge);
    }

    // Just dropped back below the threshold: log the current ridge energy
    else if ((_prev_energy >= threshold) && (_curr_energy < threshold))
    {
        if (exited_first_basin){_log_ridge_(_curr_energy);}
    }

}


RidgeEnergy::~RidgeEnergy()
{
    // Outfile closed in Base destructor
    if (energetic_threshold && !inherent_structure)
    {
        outfile = fopen(fnames.ridge_E.c_str(), "w");    
    }
    else if (!energetic_threshold && !inherent_structure)
    {
        outfile = fopen(fnames.ridge_S.c_str(), "w");    
    }
    else if (energetic_threshold && inherent_structure)
    {
        outfile = fopen(fnames.ridge_E_IS.c_str(), "w");    
    }
    else if (!energetic_threshold && inherent_structure)
    {
        outfile = fopen(fnames.ridge_S_IS.c_str(), "w");    
    }
    else
    {
        std::cout << "WARNING: data may not have been saved in ridge energy!"
        << std::endl;
    }

    fprintf(outfile, "%i %.05e %.05e %.05e %.05e %lli\n",
        int(exited_first_basin), same.mu1, var_from_S(same.S1, same.counter),
        same.current_max, same.current_min, same.counter); 
    fprintf(outfile, "%i %.05e %.05e %.05e %.05e %lli\n",
        int(exited_first_basin), diff.mu1, var_from_S(diff.S1, diff.counter),
        diff.current_max, diff.current_min, diff.counter);
}
