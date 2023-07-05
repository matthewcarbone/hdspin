#include <assert.h>

#include "simulation.h"
#include "spin.h"
#include "obs1.h"


Simulation::Simulation(const parameters::FileNames fnames,
    const parameters::SimulationParameters params) : fnames(fnames), params(params) {};

void Simulation::execute()
{
    EnergyMapping emap(params);
    SpinSystem sys(params, emap);
    parameters::StateProperties prev, curr;

    // Special case of the standard spin dynamics: if rtp.loop_dynamics == 2,
    // then the timestep is divided by rtp.N_spins.
    double waiting_time;

    // Simulation parameters
    double simulation_clock = 0.0;

    OnePointObservables obs1(fnames, params, sys);

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
        obs1.step(simulation_clock);

        if (simulation_clock > params.N_timesteps){break;}
    }
}


