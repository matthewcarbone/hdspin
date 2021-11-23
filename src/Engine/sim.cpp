/* Core local spin system algorithm.
 *
 * Matthew Carbone, Columbia University 2020
 *
 */

#include <assert.h>

#include "Utils/structures.h"
#include "Spin/gillespie.h"
#include "Spin/standard.h"
#include "Engine/sim.h"

// Observable includes
#include "Obs/energy.h"
#include "Obs/psi.h"
#include "Obs/age.h"
#include "Obs/ridge.h"


Simulation::Simulation(const FileNames fnames,
    const RuntimeParameters rtp) : fnames(fnames), rtp(rtp) {};

GillespieSimulation::GillespieSimulation(const FileNames fnames,
    const RuntimeParameters rtp) : Simulation(fnames, rtp) {};

StandardSimulation::StandardSimulation(const FileNames fnames,
    const RuntimeParameters rtp) : Simulation(fnames, rtp) {};


void GillespieSimulation::execute()
{
    GillespieSpinSystem sys(rtp);
    Vals prev, curr;
    long double waiting_time;

    Energy obs_energy(fnames);

    PsiConfig obs_psi_config(fnames, rtp, false);
    PsiConfig obs_psi_config_IS(fnames, rtp, true);

    PsiBasin obs_psi_basin_E(fnames, rtp, true, false);
    PsiBasin obs_psi_basin_S(fnames, rtp, false, false);
    PsiBasin obs_psi_basin_E_IS(fnames, rtp, true, true);
    PsiBasin obs_psi_basin_S_IS(fnames, rtp, false, true);

    AgingConfig obs_age_config(fnames);    
    AgingBasin obs_age_basin(fnames, rtp);

    // boolean arguments are inherent structure, energetic threshold, and
    // compare IS proxy.
    RidgeEnergy obs_ridge_energy_E(fnames, rtp, 0, true);
    RidgeEnergy obs_ridge_energy_S(fnames, rtp, 0, false);
    RidgeEnergy obs_ridge_energy_E_IS(fnames, rtp, 1, true);
    RidgeEnergy obs_ridge_energy_S_IS(fnames, rtp, 1, false);
    RidgeEnergy obs_ridge_energy_E_proxy_IS(fnames, rtp, 2, true);
    RidgeEnergy obs_ridge_energy_S_proxy_IS(fnames, rtp, 2, false);


    // Simulation clock is 0 before entering the while loop
    while (true)
    {
        waiting_time = sys.step_();
        simulation_clock += waiting_time;

        prev = sys.get_prev();
        curr = sys.get_curr();

        // Step observables - energy
        obs_energy.step_(simulation_clock, prev);

        // ... - psi config
        obs_psi_config.step_(waiting_time, prev, curr);
        obs_psi_config_IS.step_(waiting_time, prev, curr);

        // ... - psi basin
        obs_psi_basin_E.step_(waiting_time, prev, curr);
        obs_psi_basin_S.step_(waiting_time, prev, curr);
        obs_psi_basin_E_IS.step_(waiting_time, prev, curr);
        obs_psi_basin_S_IS.step_(waiting_time, prev, curr);

        // ... - pi config
        obs_age_config.step_(simulation_clock, sys.get_n_accept(), prev);
        
        // ... - pi basin
        obs_age_basin.step_(simulation_clock, prev, curr);

        // ... - ridge energies
        obs_ridge_energy_E.step_(prev, curr);
        obs_ridge_energy_S.step_(prev, curr);
        obs_ridge_energy_E_IS.step_(prev, curr);
        obs_ridge_energy_S_IS.step_(prev, curr);
        obs_ridge_energy_E_proxy_IS.step_(prev, curr);
        obs_ridge_energy_S_proxy_IS.step_(prev, curr);

        // Check for possible (although unlikely) overflow
        assert(simulation_clock > 0.0);

        if (simulation_clock > rtp.N_timesteps){break;}

    }
}


void StandardSimulation::execute()
{
    StandardSpinSystem sys(rtp);
    Vals prev, curr;

    // Special case of the standard spin dynamics: if rtp.loop_dynamics == 2,
    // then the timestep is divided by rtp.N_spins.
    long double waiting_time;

    // Initialize all observables
    Energy obs_energy(fnames);

    PsiConfig obs_psi_config(fnames, rtp, false);
    PsiConfig obs_psi_config_IS(fnames, rtp, true);

    PsiBasin obs_psi_basin_E(fnames, rtp, true, false);
    PsiBasin obs_psi_basin_S(fnames, rtp, false, false);
    PsiBasin obs_psi_basin_E_IS(fnames, rtp, true, true);
    PsiBasin obs_psi_basin_S_IS(fnames, rtp, false, true);

    AgingConfig obs_age_config(fnames);
    AgingBasin obs_age_basin(fnames, rtp);

    RidgeEnergy obs_ridge_energy_E(fnames, rtp, 0, true);
    RidgeEnergy obs_ridge_energy_S(fnames, rtp, 0, false);
    RidgeEnergy obs_ridge_energy_E_IS(fnames, rtp, 1, true);
    RidgeEnergy obs_ridge_energy_S_IS(fnames, rtp, 1, false);
    RidgeEnergy obs_ridge_energy_E_proxy_IS(fnames, rtp, 2, true);
    RidgeEnergy obs_ridge_energy_S_proxy_IS(fnames, rtp, 2, false);

    // Simulation clock is 0 before entering the while loop
    while (true)
    {
        
        // Standard step returns a boolean flag which is true if the new
        // proposed configuration was accepted or not.
        waiting_time = sys.step_();

        // The waiting time is always 1.0 for a standard simulation. We take
        // the convention that the "prev" structure indexes the state of the
        // spin system before the step, and that all observables are indexed
        // by the state after the step. Thus, we step the simulation_clock
        // before stepping the observables.
        simulation_clock += waiting_time;

        prev = sys.get_prev();
        curr = sys.get_curr();

        // Step observables - energy
        obs_energy.step_(simulation_clock, prev);

        // ... - psi config
        obs_psi_config.step_(waiting_time, prev, curr);
        obs_psi_config_IS.step_(waiting_time, prev, curr);

        // ... - psi basin
        obs_psi_basin_E.step_(waiting_time, prev, curr);
        obs_psi_basin_S.step_(waiting_time, prev, curr);
        obs_psi_basin_E_IS.step_(waiting_time, prev, curr);
        obs_psi_basin_S_IS.step_(waiting_time, prev, curr);

        // ... - pi config
        obs_age_config.step_(simulation_clock, sys.get_n_accept(), prev);
        
        // ... - pi basin
        obs_age_basin.step_(simulation_clock, prev, curr);

        // ... - ridge energies
        obs_ridge_energy_E.step_(prev, curr);
        obs_ridge_energy_S.step_(prev, curr);
        obs_ridge_energy_E_IS.step_(prev, curr);
        obs_ridge_energy_S_IS.step_(prev, curr);
        obs_ridge_energy_E_proxy_IS.step_(prev, curr);
        obs_ridge_energy_S_proxy_IS.step_(prev, curr);

        if (simulation_clock > rtp.N_timesteps){break;}

    }
}
