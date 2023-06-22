#include <assert.h>

#include "simulation.h"
#include "spin.h"
#include "obs_energy.h"


Simulation::Simulation(const parameters::FileNames fnames,
    const parameters::SimulationParameters params) : fnames(fnames), params(params) {};

// GillespieSimulation::GillespieSimulation(const parameters::FileNames fnames,
//     const parameters::SimulationParameters params) : Simulation(fnames, params) {};

StandardSimulation::StandardSimulation(const parameters::FileNames fnames,
    const parameters::SimulationParameters params) : Simulation(fnames, params) {};


// void GillespieSimulation::execute()
// {

//     EnergyMapping emap(params);
//     GillespieSpinSystem sys(params, emap);

//     Vals prev, curr;
//     long double waiting_time;

//     Energy obs_energy(fnames, rtp);
//     EnergyAvgNeighbors obs_energy_avg_neighbors(fnames, rtp);
//     EnergyInherentStructure obs_energy_IS(fnames, rtp);
//     Configs obs_config(fnames, rtp);

//     PsiConfig obs_psi_config(fnames, rtp);
//     PsiConfigInherentStructure obs_psi_config_IS(fnames, rtp);

//     PsiBasinThreshold obs_psi_basin_E(fnames, rtp);
//     PsiBasinAttractor obs_psi_basin_S(fnames, rtp);
//     PsiBasinThresholdInherentStructure obs_psi_basin_E_IS(fnames, rtp);
//     PsiBasinAttractorInherentStructure obs_psi_basin_S_IS(fnames, rtp);

//     AgingConfig obs_age_config(fnames, rtp);
//     AgingConfigInherentStructure obs_age_config_IS(fnames, rtp);
//     AgingBasinThreshold obs_age_basin_E(fnames, rtp);
//     AgingBasinThresholdInherentStructure obs_age_basin_E_IS(fnames, rtp);
//     AgingBasinAttractor obs_age_basin_S(fnames, rtp);
//     AgingBasinAttractorInherentStructure obs_age_basin_S_IS(fnames, rtp);

//     // boolean arguments are inherent structure, energetic threshold, and
//     // compare IS proxy.
//     RidgeEnergy obs_ridge_E(fnames, rtp);
//     RidgeAttractor obs_ridge_S(fnames, rtp);

//     // Simulation clock is 0 before entering the while loop
//     while (true)
//     {
//         waiting_time = sys.step();

//         simulation_clock += waiting_time;

//         prev = sys.get_prev();
//         curr = sys.get_curr();

//         // Step observables - energy
//         step_all_observables(prev, curr, waiting_time, simulation_clock, obs_energy);
//         obs_energy_avg_neighbors.step(simulation_clock, sys);
//         obs_energy_IS.step(simulation_clock, prev);
//         obs_config.step(simulation_clock, prev);

//         // ... - psi config
//         obs_psi_config.step(waiting_time, prev, curr);
//         obs_psi_config_IS.step(waiting_time, prev, curr);

//         // ... - psi basin
//         obs_psi_basin_E.step(waiting_time, prev, curr);
//         obs_psi_basin_S.step(waiting_time, prev, curr);
//         obs_psi_basin_E_IS.step(waiting_time, prev, curr);
//         obs_psi_basin_S_IS.step(waiting_time, prev, curr);

//         // ... - pi config
//         obs_age_config.step(simulation_clock, prev);
//         obs_age_config_IS.step(simulation_clock, prev);
        
//         // ... - pi basin
//         obs_age_basin_E.step(simulation_clock, prev, curr);
//         obs_age_basin_S.step(simulation_clock, prev, curr);
//         obs_age_basin_E_IS.step(simulation_clock, prev, curr);
//         obs_age_basin_S_IS.step(simulation_clock, prev, curr);

//         // ... - ridge energies
//         obs_ridge_E.step(prev, curr, waiting_time, simulation_clock);
//         obs_ridge_S.step(prev, curr, waiting_time, simulation_clock);

//         // Check for possible (although unlikely) overflow
//         assert(simulation_clock > 0.0);

//         if (simulation_clock > rtp.N_timesteps)
//         {
//             break;
//         }
//     }
// }


void StandardSimulation::execute()
{
    EnergyMapping emap(params);
    StandardSpinSystem sys(params, emap);
    parameters::StateProperties prev, curr;

    // Special case of the standard spin dynamics: if rtp.loop_dynamics == 2,
    // then the timestep is divided by rtp.N_spins.
    long double waiting_time;

    // Simulation parameters
    long double simulation_clock = 0.0;

    Energy obs_energy(fnames, params);

    // Simulation clock is 0 before entering the while loop
    while (true)
    {
        
        // Standard step returns a boolean flag which is true if the new
        // proposed configuration was accepted or not.
        waiting_time = sys.step();

        // The waiting time is always 1.0 for a standard simulation. We take
        // the convention that the "prev" structure indexes the state of the
        // spin system before the step, and that all observables are indexed
        // by the state after the step. Thus, we step the simulation_clock
        // before stepping the observables.
        simulation_clock += waiting_time;

        prev = sys.get_previous_state();
        curr = sys.get_current_state();

        // Step observables
        obs_energy.step(simulation_clock, prev.energy);

        if (simulation_clock > params.N_timesteps - 1){break;}
    }

    sys.summarize();
}
