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
#include "obs2.h"
#include "CLI11/CLI11.hpp"


void step_all_observables_(const double waiting_time, const double simulation_clock, OnePointObservables& obs1, RidgeE& ridgeE, RidgeS& ridgeS, PsiConfig& psiConfig)
{
    obs1.step(waiting_time, simulation_clock);
    ridgeE.step(waiting_time, simulation_clock);
    ridgeS.step(waiting_time, simulation_clock);
    psiConfig.step(waiting_time);
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
    PsiConfig psiConfig(fnames, params, sys);

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
        // before stepping the observables. Note that the waiting time can
        // vary for the Gillespie dynamics.
        simulation_clock += waiting_time;

        step_all_observables_(
            waiting_time,
            simulation_clock,
            obs1,
            ridgeE,
            ridgeS,
            psiConfig
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
    for (unsigned int step=0; step<int(1e7); step++)
    {
        simulation_clock += sys.step();
        if (simulation_clock > p.N_timesteps){break;}
    }
    return time_utils::get_time_delta(t_start) / simulation_clock;
}

std::string determine_dynamics_automatically(const parameters::SimulationParameters params, const unsigned int mpi_world_size, const unsigned int mpi_rank, MPI_Comm mpi_comm)
{
    double standard_time = 0.0;
    double gillespie_time = 0.0;
    double standard_std = 0.0;
    double gillespie_std = 0.0;
    parameters::SimulationParameters p = params;
    unsigned int result_int = 2;

    // Run both on rank 1 if we only have a single process
    if (mpi_world_size == 1)
    {
        standard_time = get_sim_time(p, "standard");
        gillespie_time = get_sim_time(p, "gillespie");
    }

    // Otherwise, we actually want to divide up the work a bit
    // Let all even ranks (including 0) run the standard simulation
    // and all odd ranks run the Gillespie simulation
    // The results can then be averaged at the end by ranks 0 and 1.
    else
    {
        double times[mpi_world_size];
        if (mpi_rank % 2 == 0)
        {
            times[mpi_rank] = get_sim_time(p, "standard");
            // MPI_Recv(&gillespie_time, 1, MPI_DOUBLE, 1, 0, mpi_comm, MPI_STATUS_IGNORE);
        }
        else if (mpi_rank % 2 != 0)
        {
            times[mpi_rank] = get_sim_time(p, "gillespie");
            // MPI_Send(&gillespie_time, 1, MPI_DOUBLE, 0, 0, mpi_comm);
        }
        MPI_Barrier(mpi_comm);

        // Now, we send everything to rank 0
        if (mpi_rank != 0)
        {
            MPI_Send(&times[mpi_rank], 1, MPI_DOUBLE, 0, 0, mpi_comm);
        }
        else
        {
            for (unsigned int rank=1; rank<mpi_world_size; rank++)
            {
                MPI_Recv(&times[rank], 1, MPI_DOUBLE, rank, 0, mpi_comm, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(mpi_comm);

        // And let rank 0 deal with the rest
        if (mpi_rank == 0)
        {
            std::vector<double> standard_times;
            std::vector<double> gillespie_times;
            for (unsigned int ii=0; ii<mpi_world_size; ii++)
            {
                if (ii % 2 == 0){standard_times.push_back(times[ii]);}
                else{gillespie_times.push_back(times[ii]);}
            }

            // Calculate the mean and standard deviation
            standard_time = mean_vector(standard_times);
            standard_std = sqrt(variance_vector(standard_times));
            gillespie_time = mean_vector(gillespie_times);
            gillespie_std = sqrt(variance_vector(gillespie_times));
        }
    }

    if (mpi_rank == 0)
    {
        printf("Gillespie vs. Standard dynamics:\n\t%.02e +/- %.02e vs. %.02e +/- %.02e wall/sim\n", gillespie_time, gillespie_std, standard_time, standard_std);
        if (gillespie_time < standard_time)
        {
            printf("\tRunning Gillespie dynamics, faster by factor of %.01f\n", standard_time / gillespie_time);
            result_int = 1;
        }
        else
        {
            printf("\tRunning standard dynamics, faster by factor of %.01f\n", gillespie_time / standard_time);   
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

    // Define the defaults
    parameters::SimulationParameters p;

    // PARSER ----------------------------------------------------------------
    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    CLI::App app{
        "hdspin is an application for simulating binary spin systems"
    };

    app.add_option(
        "-N, --N_spins", p.N_spins,
        "Number of binary spins to use in the simulation. Must be bounded "
        "between 1 and the PRECISON option set at compile time."
    )->check(CLI::Range(1, PRECISON))->required();

    app.add_option(
        "-l, --landscape", p.landscape,
        "Choice of landscape, either the exponential random energy model "
        "(EREM) or the Gaussian random energy model (GREM)."
    )->check(CLI::IsMember({"EREM", "GREM"}))->required();

    app.add_option(
        "-b, --beta", p.beta,
        "Inverse temperature."
    )->check(CLI::PositiveNumber)->required();

    app.add_option(
        "-t, --log10_N_timesteps", p.log10_N_timesteps,
        "The number of timesteps to run on a log10 scale. For example, if "
        "'t=7', then a simulation of length 10^7 will be run."
    )->check(CLI::PositiveNumber)->required();

    app.add_option(
        "-m, --memory", p.memory,
        "The size of the memory cache. If -1, will attempt to use a cache "
        "size equal to 2^N_spins, which might cause memory issues. The "
        "default is 2^25."
    )->check(CLI::PositiveNumber|CLI::IsMember({-1}));

    app.add_option(
        "-d, --dynamics", p.dynamics,
        "The type of dynamics to run. Defaults to 'auto'. Standard dynamics "
        "attempts to flip one spin every timestep, with the Metropolis "
        "acceptance/rejection criterion. Gillespie dynamics calculates all "
        "exit rates at once and flips a spin every iteration of the "
        "algorithm, but with a waiting time not necessarily equal to 1. The "
        "auto selection runs quick simulations of both types to see which "
        "is faster, and selects that one."
    )->check(CLI::IsMember({"standard", "gillespie", "auto"}));

    app.add_option(
        "-n, --n_tracers_per_MPI_rank", p.n_tracers_per_MPI_rank,
        "The number of simulations per MPI rank to run. Defaults to 10."
    )->check(CLI::PositiveNumber);

    app.add_option(
        "--seed", p.seed,
        "Seeds for the random number generators for reproducible runs. Leave "
        "unset for un-reproducible, random seeds."
    )->check(CLI::PositiveNumber);

    app.add_flag("--skip_IS{false}", p.calculate_inherent_structure_observables,
        "Providing the --skip_IS flag informs hdspin to skip all calculations of the inherent structure. This will speed up the simulation and use ~log N less memory in the cache.");

    CLI11_PARSE(app, argc, argv);
    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    // PARSER ----------------------------------------------------------------

    update_parameters_(&p);

    MPI_Barrier(MPI_COMM_WORLD);

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

    // Get the information for this MPI rank
    const unsigned int n_tracers_per_MPI_rank = p.n_tracers_per_MPI_rank;
    const unsigned int resume_at = 0;
    const unsigned int start = resume_at + MPI_RANK * n_tracers_per_MPI_rank;
    const unsigned int end = resume_at + (MPI_RANK + 1) * n_tracers_per_MPI_rank;

    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_RANK == 0)
    {
        // parameters::log_json(inp);
        parameters::log_parameters(p);
        make_directories();
        grids::make_energy_grid_logspace(p.log10_N_timesteps, p.grid_size);
        grids::make_pi_grids(p.log10_N_timesteps, p.dw, p.grid_size);
        json j = parameters::parameters_to_json(p);
        std::ofstream o("config.json");
        o << std::setw(4) << j << std::endl;
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
