#include <algorithm>    // std::min

#include "Obs/base.h"
#include "Obs/rolling.h"
#include "Utils/utils.h"
#include "Utils/structures.h"

Rolling::Rolling(const FileNames fnames, const RuntimeParameters rtp) :
    Base(fnames), rtp(rtp) {}


void Rolling::_log_ridge_E(const Vals curr)
{
    // Check the condition for whether or not the last energy before the
    // ridge is the same as the current energy.

    if (curr.energy == E_last_energy)
    {
        E_e_same.mu1 = iterative_mean(E_e_same.mu0,
            E_current_ridge, E_e_same.counter);
        E_e_same.S1 = iterative_S(E_e_same.mu0, E_e_same.mu1, E_current_ridge,
            E_e_same.S0);
        E_e_same.mu0 = E_e_same.mu1;
        E_e_same.S0 = E_e_same.S1;
        E_e_same.current_min = std::min(E_e_same.current_min, E_current_ridge);
        E_e_same.current_max = std::max(E_e_same.current_max, E_current_ridge);
        E_e_same.counter += 1;
    }
    else
    {
        E_e_diff.mu1 = iterative_mean(E_e_diff.mu0,
            E_current_ridge, E_e_diff.counter);
        E_e_diff.S1 = iterative_S(E_e_diff.mu0, E_e_diff.mu1, E_current_ridge,
            E_e_diff.S0);
        E_e_diff.mu0 = E_e_diff.mu1;
        E_e_diff.S0 = E_e_diff.S1;
        E_e_diff.current_min = std::min(E_e_diff.current_min, E_current_ridge);
        E_e_diff.current_max = std::max(E_e_diff.current_max, E_current_ridge);
        E_e_diff.counter += 1;
    }
}



void Rolling::_log_ridge_E_IS(const Vals curr)
{
    // Check the condition for whether or not the last energy before the
    // ridge is the same as the current energy.

    if (curr.energy_IS == E_IS_last_energy)
    {
        E_IS_e_same.mu1 = iterative_mean(E_IS_e_same.mu0,
            E_IS_current_ridge, E_IS_e_same.counter);
        E_IS_e_same.S1 = iterative_S(E_IS_e_same.mu0, E_IS_e_same.mu1, 
            E_IS_current_ridge, E_IS_e_same.S0);
        E_IS_e_same.mu0 = E_IS_e_same.mu1;
        E_IS_e_same.S0 = E_IS_e_same.S1;
        E_IS_e_same.current_min = std::min(E_IS_e_same.current_min,
            E_IS_current_ridge);
        E_IS_e_same.current_max = std::max(E_IS_e_same.current_max,
            E_IS_current_ridge);
        E_IS_e_same.counter += 1;
    }
    else
    {
        E_IS_e_diff.mu1 = iterative_mean(E_IS_e_diff.mu0,
            E_IS_current_ridge, E_IS_e_diff.counter);
        E_IS_e_diff.S1 = iterative_S(E_IS_e_diff.mu0,
            E_IS_e_diff.mu1, E_IS_current_ridge,
            E_IS_e_diff.S0);
        E_IS_e_diff.mu0 = E_IS_e_diff.mu1;
        E_IS_e_diff.S0 = E_IS_e_diff.S1;
        E_IS_e_diff.current_min = std::min(E_IS_e_diff.current_min,
            E_IS_current_ridge);
        E_IS_e_diff.current_max = std::max(E_IS_e_diff.current_max,
            E_IS_current_ridge);
        E_IS_e_diff.counter += 1;
    }
}


void Rolling::_log_ridge_S(const Vals curr)
{
    // Check the condition for whether or not the last energy before the
    // ridge is the same as the current energy.

    if (curr.energy == S_last_energy)
    {
        S_e_same.mu1 = iterative_mean(S_e_same.mu0,
            S_current_ridge, S_e_same.counter);
        S_e_same.S1 = iterative_S(S_e_same.mu0, S_e_same.mu1, S_current_ridge,
            S_e_same.S0);
        S_e_same.mu0 = S_e_same.mu1;
        S_e_same.S0 = S_e_same.S1;
        S_e_same.current_min = std::min(S_e_same.current_min, S_current_ridge);
        S_e_same.current_max = std::max(S_e_same.current_max, S_current_ridge);
        S_e_same.counter += 1;
    }
    else
    {
        S_e_diff.mu1 = iterative_mean(S_e_diff.mu0,
            S_current_ridge, S_e_diff.counter);
        S_e_diff.S1 = iterative_S(S_e_diff.mu0, S_e_diff.mu1, S_current_ridge,
            S_e_diff.S0);
        S_e_diff.mu0 = S_e_diff.mu1;
        S_e_diff.S0 = S_e_diff.S1;
        S_e_diff.current_min = std::min(S_e_diff.current_min, S_current_ridge);
        S_e_diff.current_max = std::max(S_e_diff.current_max, S_current_ridge);
        S_e_diff.counter += 1;
    }
}


void Rolling::_log_ridge_S_IS(const Vals curr)
{
    // Check the condition for whether or not the last energy before the
    // ridge is the same as the current energy.

    if (curr.energy_IS == S_IS_last_energy)
    {
        S_IS_e_same.mu1 = iterative_mean(S_IS_e_same.mu0,
            S_IS_current_ridge, S_IS_e_same.counter);
        S_IS_e_same.S1 = iterative_S(S_IS_e_same.mu0, S_IS_e_same.mu1, 
            S_IS_current_ridge, S_IS_e_same.S0);
        S_IS_e_same.mu0 = S_IS_e_same.mu1;
        S_IS_e_same.S0 = S_IS_e_same.S1;
        S_IS_e_same.current_min = std::min(S_IS_e_same.current_min,
            S_IS_current_ridge);
        S_IS_e_same.current_max = std::max(S_IS_e_same.current_max,
            S_IS_current_ridge);
        S_IS_e_same.counter += 1;
    }
    else
    {
        S_IS_e_diff.mu1 = iterative_mean(S_IS_e_diff.mu0,
            S_IS_current_ridge, S_IS_e_diff.counter);
        S_IS_e_diff.S1 = iterative_S(S_IS_e_diff.mu0,
            S_IS_e_diff.mu1, S_IS_current_ridge,
            S_IS_e_diff.S0);
        S_IS_e_diff.mu0 = S_IS_e_diff.mu1;
        S_IS_e_diff.S0 = S_IS_e_diff.S1;
        S_IS_e_diff.current_min = std::min(S_IS_e_diff.current_min,
            S_IS_current_ridge);
        S_IS_e_diff.current_max = std::max(S_IS_e_diff.current_max,
            S_IS_current_ridge);
        S_IS_e_diff.counter += 1;
    }
}


void Rolling::step_(const Vals prev, const Vals curr)
{

    // Handle the energetic barrier
    if ((prev.energy < rtp.energetic_threshold) & (curr.energy >=
        rtp.energetic_threshold))
    {
        // Just exited a basin, track the previous energy
        E_last_energy = prev.energy;
        E_current_ridge = curr.energy;
    }
    else if ((prev.energy >= rtp.energetic_threshold) & (curr.energy >=
        rtp.energetic_threshold))
    {
        // Still above
        E_current_ridge = std::max(curr.energy, E_current_ridge);
    }
    else if ((prev.energy >= rtp.energetic_threshold) & (curr.energy <
        rtp.energetic_threshold))
    {
        // Just dropped back below the threshold
        // Want to log the current ridge energy
        _log_ridge_E(curr);
    }

    // Handle the entropic attractor
    if ((prev.energy < rtp.entropic_attractor) & (curr.energy >=
        rtp.entropic_attractor))
    {
        // Just exited a basin, track the previous energy
        S_last_energy = prev.energy;
        S_current_ridge = curr.energy;
    }
    else if ((prev.energy >= rtp.entropic_attractor) & (curr.energy >=
        rtp.entropic_attractor))
    {
        // Still above
        S_current_ridge = std::max(curr.energy, S_current_ridge);
    }
    else if ((prev.energy >= rtp.entropic_attractor) & (curr.energy <
        rtp.entropic_attractor))
    {
        // Just dropped back below the threshold
        // Want to log the current ridge energy
        _log_ridge_S(curr);
    }

    // Same thing for the inherent structure -------


    // Handle the energetic barrier
    if ((prev.energy_IS < rtp.energetic_threshold) & (curr.energy_IS >=
        rtp.energetic_threshold))
    {
        // Just exited a basin, track the previous energy
        E_IS_last_energy = prev.energy_IS;
        E_IS_current_ridge = curr.energy_IS;
    }
    else if ((prev.energy_IS >= rtp.energetic_threshold) & (curr.energy_IS >=
        rtp.energetic_threshold))
    {
        // Still above
        E_IS_current_ridge = std::max(curr.energy_IS, E_IS_current_ridge);
    }
    else if ((prev.energy_IS >= rtp.energetic_threshold) & (curr.energy_IS <
        rtp.energetic_threshold))
    {
        // Just dropped back below the threshold
        // Want to log the current ridge energy
        _log_ridge_E_IS(curr);
    }

    // Handle the entropic attractor
    if ((prev.energy_IS < rtp.entropic_attractor) & (curr.energy_IS >=
        rtp.entropic_attractor))
    {
        // Just exited a basin, track the previous energy
        S_IS_last_energy = prev.energy_IS;
        S_IS_current_ridge = curr.energy_IS;
    }
    else if ((prev.energy_IS >= rtp.entropic_attractor) & (curr.energy_IS >=
        rtp.entropic_attractor))
    {
        // Still above
        S_IS_current_ridge = std::max(curr.energy_IS, S_IS_current_ridge);
    }
    else if ((prev.energy_IS >= rtp.entropic_attractor) & (curr.energy_IS <
        rtp.entropic_attractor))
    {
        // Just dropped back below the threshold
        // Want to log the current ridge energy
        _log_ridge_S_IS(curr);
    }
}


Rolling::~Rolling()
{
    // Outfile closed in Base destructor
    outfile = fopen(fnames.rolling.c_str(), "w");

    // Dump the results to disk for the standard trajectory
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        E_e_same.mu1, var_from_S(E_e_same.S1, E_e_same.counter),
        E_e_same.current_max, E_e_same.current_min, E_e_same.counter); 
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        E_e_diff.mu1, var_from_S(E_e_diff.S1, E_e_diff.counter),
        E_e_diff.current_max, E_e_diff.current_min, E_e_diff.counter);
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        S_e_same.mu1, var_from_S(S_e_same.S1, S_e_same.counter),
        S_e_same.current_max, S_e_same.current_min, S_e_same.counter);
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        S_e_diff.mu1, var_from_S(S_e_diff.S1, S_e_diff.counter),
        S_e_diff.current_max, S_e_diff.current_min, S_e_diff.counter);

    // Dump the inherent structure trajectory results
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        E_IS_e_same.mu1, var_from_S(E_IS_e_same.S1, E_IS_e_same.counter),
        E_IS_e_same.current_max, E_IS_e_same.current_min, E_IS_e_same.counter);
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        E_IS_e_diff.mu1, var_from_S(E_IS_e_diff.S1, E_IS_e_diff.counter),
        E_IS_e_diff.current_max, E_IS_e_diff.current_min, E_IS_e_diff.counter);
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        S_IS_e_same.mu1, var_from_S(S_IS_e_same.S1, S_IS_e_same.counter),
        S_IS_e_same.current_max, S_IS_e_same.current_min, S_IS_e_same.counter);
    fprintf(outfile, "%.05Le %.05Le %.05Le %.05Le %lli\n",
        S_IS_e_diff.mu1, var_from_S(S_IS_e_diff.S1, S_IS_e_diff.counter),
        S_IS_e_diff.current_max, S_IS_e_diff.current_min, S_IS_e_diff.counter);
}
