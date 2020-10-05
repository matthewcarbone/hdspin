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
    PsiConfig obs_psi_config(fnames, rtp);
    AgingConfig obs_age_config(fnames);
    PsiBasin obs_psi_basin(fnames, rtp);
    AgingBasin obs_age_basin(fnames, rtp);

    // Simulation clock is 0 before entering the while loop
    while (true)
    {
        waiting_time = sys.step_();
        simulation_clock += waiting_time;

        prev = sys.get_prev();
        curr = sys.get_curr();

        // Step observables
        obs_energy.step_(simulation_clock, prev);
        obs_psi_config.step_(waiting_time, prev, curr);
        obs_age_config.step_(simulation_clock, sys.get_n_accept(), prev);
        obs_psi_basin.step_(waiting_time, prev, curr);
        obs_age_basin.step_(simulation_clock, prev, curr);


        /*
        //         -------------------------------------------------
        //         ---------------- STEP TRACKERS ------------------
        //         -------------------------------------------------

        // Append trackers. Note that Gillespie dynamics are different than
        // standard in the order in which we update the grids, so the grids are
        // actually stepped before the sys/inh objects are updated.  
        energy_grid.step(current_time, sys, inh);
        aging_config_grid.step(current_time, n_accept, sys.x, inh.x);
        psi_config_counter.step(waiting_time, false);  // Step standard

        // This is a tricky update for the inherent structure, since it
        // will have a different waiting time than the normal
        // configuration, as it may not change even though the normal
        // configuration does.
        if (inh.x == inh.x_prev){inh.waiting_time += waiting_time;}
        else
        {
            // Step the inherent structure psi config counter
            psi_config_counter.step(inh.waiting_time, true);  // Step IS
            inh.waiting_time = 0.0;
        }

        psi_basin_counter.step(&sys, &inh);

        //         -------------------------------------------------
        //         -------------- DONE STEP TRACKERS ---------------
        //         -------------------------------------------------

        */
        // --------------------------------------------------------------------
        // ----------------------- ENGINE FINISH ------------------------------
        // --------------------------------------------------------------------

        // Check for possible (although unlikely) overflow
        assert(simulation_clock > 0.0);

        if (simulation_clock > rtp.N_timesteps){break;}

    }

    // Close the outfiles and write to disk when not doing so dynamically
    /*
    energy_grid.close_outfile();
    psi_config_counter.write_to_disk(fnames.psi_config);
    psi_basin_counter.write_to_disk(fnames.psi_basin);
    aging_config_grid.close_outfile();
    */
}


void StandardSimulation::execute()
{
    StandardSpinSystem sys(rtp);
    Vals prev, curr;
    bool accepted;
    const long double waiting_time = 1.0;

    Energy obs_energy(fnames);
    PsiConfig obs_psi_config(fnames, rtp);
    AgingConfig obs_age_config(fnames);
    PsiBasin obs_psi_basin(fnames, rtp);
    AgingBasin obs_age_basin(fnames, rtp);

    // Simulation clock is 0 before entering the while loop
    while (true)
    {
        
        // Standard step returns a boolean flag which is true if the new
        // proposed configuration was accepted or not.
        accepted = sys.step_();

        // The waiting time is always 1.0 for a standard simulation. We take
        // the convention that the "prev" structure indexes the state of the
        // spin system before the step, and that all observables are indexed
        // by the state after the step. Thus, we step the simulation_clock
        // before stepping the observables.
        simulation_clock += waiting_time;

        prev = sys.get_prev();
        curr = sys.get_curr();

        // Step observables
        obs_energy.step_(simulation_clock, prev);
        obs_psi_config.step_(waiting_time, prev, curr);
        obs_age_config.step_(simulation_clock, sys.get_n_accept(), prev);
        obs_psi_basin.step_(waiting_time, prev, curr);
        obs_age_basin.step_(simulation_clock, prev, curr);


        // --------------------------------------------------------------------
        // ----------------------- ENGINE FINISH ------------------------------
        // --------------------------------------------------------------------

        //         -------------------------------------------------
        //         ------------- STEP [other] TRACKERS -------------
        //         -------------------------------------------------

        /*
        energy_grid.step(timestep, sys, inh);
        aging_config_grid.step(timestep, n_accepted, sys.x, inh.x);
        update_basin_information(&sys, params, 1.0);
        update_basin_information(&inh, params, 1.0);
        psi_basin_counter.step(&sys, &inh);
        */

        //         -------------------------------------------------
        //         -------------- DONE STEP TRACKERS ---------------
        //         -------------------------------------------------

        if (simulation_clock > rtp.N_timesteps){break;}

    }

    /*
    // Close the outfiles and write to disk when not doing so dynamically
    
    psi_config_counter.write_to_disk(fnames.psi_config);
    psi_basin_counter.write_to_disk(fnames.psi_basin);
    aging_config_grid.close_outfile();

    delete[] energy_arr;
    delete[] inherent_structure_mapping;
    */
}
