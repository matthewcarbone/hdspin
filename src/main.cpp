#include <fstream>      // std::ofstream
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <sstream>
#include <assert.h>
#include <vector>
#include <set>
#include <mpi.h>

#include "utils.h"
#include "spin.h"
#include "obs1.h"


void step_all_observables_(const double waiting_time, const double simulation_clock, OnePointObservables& obs1, RidgeE& ridgeE, RidgeS& ridgeS)
{
    obs1.step(waiting_time, simulation_clock);
    ridgeE.step(waiting_time, simulation_clock);
    ridgeS.step(waiting_time, simulation_clock);
}

void execute(const parameters::FileNames fnames,
    const parameters::SimulationParameters params)
{
    EnergyMapping emap(params);
    SpinSystem sys(params, emap);

    // Special case of the standard spin dynamics: if rtp.loop_dynamics == 2,
    // then the timestep is divided by rtp.N_spins.
    double waiting_time;

    // Simulation parameters
    double simulation_clock = 0.0;

    RidgeE ridgeE(fnames, params, sys);
    RidgeS ridgeS(fnames, params, sys);
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

        step_all_observables_(
            waiting_time,
            simulation_clock,
            obs1,
            ridgeE,
            ridgeS
        );

        if (simulation_clock > params.N_timesteps){break;}
    }
}

double get_sim_time(parameters::SimulationParameters p, const std::string dynamics)
{
    double simulation_clock = 0.0;
    p.dynamics = dynamics;
    EnergyMapping emap(p);
    SpinSystem sys(p, emap);
    auto t_start = std::chrono::high_resolution_clock::now();
    while (true)
    {
        simulation_clock += sys.step();
        if (simulation_clock > p.N_timesteps){break;}
    }
    return time_utils::get_time_delta(t_start);
}

std::string determine_dynamics_automatically(const parameters::SimulationParameters params, const unsigned int mpi_world_size, const unsigned int mpi_rank, MPI_Comm mpi_comm)
{
    double standard_time, gillespie_time, simulation_clock;
    parameters::SimulationParameters p = params;
    unsigned int result_int = 2;

    // Run both on rank 1
    if (mpi_world_size == 1)
    {
        standard_time = get_sim_time(p, "standard");
        gillespie_time = get_sim_time(p, "gillespie");
    }
    else
    {
        if (mpi_rank == 0)
        {
            standard_time = get_sim_time(p, "standard");
            MPI_Recv(&gillespie_time, 1, MPI_DOUBLE, 1, 0, mpi_comm, MPI_STATUS_IGNORE);
        }
        else if (mpi_rank == 1)
        {
            gillespie_time = get_sim_time(p, "gillespie");
            MPI_Send(&gillespie_time, 1, MPI_DOUBLE, 0, 0, mpi_comm);
        }
        MPI_Barrier(mpi_comm);
    }

    if (mpi_rank == 0)
    {
        printf("Gillespie vs. Standard dynamics: %.02e vs. %.02e s/simulation\n", gillespie_time, standard_time);
        if (gillespie_time < standard_time)
        {
            printf("Running Gillespie dynamics, faster by factor of %.01f\n", standard_time / gillespie_time);
            result_int = 1;
        }
        else
        {
            printf("Running standard dynamics, faster by factor of %.01f\n", gillespie_time / standard_time);   
            result_int = 0;
        }
    }

    MPI_Bcast(&result_int, 1, MPI_INT, 0, mpi_comm);
    MPI_Barrier(mpi_comm);
    
    if (result_int == 1){return "gillespie";}
    else if (result_int == 0){return "standard";}
    else{MPI_Abort(mpi_comm, 1); return "error";}
}

int main(int argc, char *argv[])
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int MPI_WORLD_SIZE;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_WORLD_SIZE);

    // Get the rank of the process
    int MPI_RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Ready: processor %s, rank %d/%d\n", processor_name, MPI_RANK, MPI_WORLD_SIZE);

    // Quick barrier to make sure the printing works out cleanly
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // Load the parameters
    std::string params_name;
    if (argc == 1){params_name = "config.json";}
    else{params_name = argv[1];}
    if (MPI_RANK == 0){printf("Loading config from %s", params_name.c_str());}

    std::ifstream ifs(params_name);
    json inp = json::parse(ifs);
    parameters::SimulationParameters p = parameters::get_parameters(inp);

    MPI_Barrier(MPI_COMM_WORLD);

    // Get the information for this MPI rank
    const unsigned int n_tracers_per_MPI_rank = p.n_tracers_per_MPI_rank;
    const unsigned int resume_at = 0;
    const unsigned int start = resume_at + MPI_RANK * n_tracers_per_MPI_rank;
    const unsigned int end = resume_at + (MPI_RANK + 1) * n_tracers_per_MPI_rank;

    // So everything prints cleanly
    printf("RANK %i/%i ID's %i -> %i\n", MPI_RANK, MPI_WORLD_SIZE, start, end);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_RANK == 0)
    {
        parameters::log_json(inp);
        parameters::log_parameters(p);
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_RANK == 0)
    {
        make_directories();
        grids::make_energy_grid_logspace(p.log10_N_timesteps, p.grid_size);
        grids::make_pi_grids(p.log10_N_timesteps, p.dw, p.grid_size);
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // Define some helpers to be used to track progress.
    const unsigned int total_steps = end - start;
    unsigned int step_size = total_steps / 10; // Print at 10 percent steps
    if (step_size == 0){step_size = 1;}
    unsigned int loop_count = 0;

    auto global_start = std::chrono::high_resolution_clock::now();

    const unsigned int starting_seed = p.seed;

    // If p.dynamics == "auto", we run a quick check to see which
    // simulation is faster
    if (p.dynamics == "auto")
    {
        p.dynamics = determine_dynamics_automatically(p, MPI_WORLD_SIZE, MPI_RANK, MPI_COMM_WORLD);
    }

    for(int ii=start; ii<end; ii++)
    {

        auto t_start = std::chrono::high_resolution_clock::now();

        const parameters::FileNames fnames = parameters::get_filenames(ii);

        // Change the seed based on the MPI rank, very important for seeded runs!
        // This will be ignored later if p.use_manual_seed is false
        p.seed = starting_seed + ii + MPI_RANK * n_tracers_per_MPI_rank;

        // Run dynamics START -------------------------------------------------
        execute(fnames, p);
        // Run dynamics END ---------------------------------------------------

        const double duration = time_utils::get_time_delta(t_start);

        loop_count++;

        if (MPI_RANK == 0)
        {
            if (loop_count % step_size == 0 | loop_count == 1)
            {
                const std::string dt_string = time_utils::get_datetime();
                const double global_duration = time_utils::get_time_delta(global_start);
                
                printf(
                    "%s ~ %s done in %.01f s (%i/%i) total elapsed %.01f s\n", dt_string.c_str(), fnames.ii_str.c_str(), duration, loop_count, total_steps, global_duration
                );
                fflush(stdout);
            }
        }
    }

    MPI_Finalize();
}
