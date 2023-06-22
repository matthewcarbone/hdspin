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
#include "simulation.h"


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
    printf("Ready: processor %s, rank %d/%d\n", processor_name, MPI_RANK,
        MPI_WORLD_SIZE);

    // Quick barrier to make sure the printing works out cleanly
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // parameters::SimulationParameters p = parameters::get_parameters();

    //////-------///////-------///////-------///////-------///////-------//////
    //////-------///////-------///////-------///////-------///////-------//////
    //////-------///////-------///////-------///////-------///////-------//////
    // Define some testing simulation parameters
    parameters::SimulationParameters p;
    p.log10_N_timesteps = 6;
    p.N_timesteps = ipow(10, int(p.log10_N_timesteps));
    p.N_spins = 10;
    p.landscape = "EREM";
    p.beta = 2.4;
    p.beta_critical = 1.0;
    p.dynamics = "standard";
    p.memory = 15;
    p.n_tracers_per_MPI_rank = 1;
    //////-------///////-------///////-------///////-------///////-------//////
    //////-------///////-------///////-------///////-------///////-------//////
    //////-------///////-------///////-------///////-------///////-------//////

    // Get the information for this MPI rank
    const int n_tracers_per_MPI_rank = p.n_tracers_per_MPI_rank;
    const int resume_at = 0;
    const int start = resume_at + MPI_RANK * n_tracers_per_MPI_rank;
    const int end = resume_at + (MPI_RANK + 1) * n_tracers_per_MPI_rank;

    // So everything prints cleanly
    printf("RANK %i/%i ID's %i -> %i\n", MPI_RANK, MPI_WORLD_SIZE, start, end);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_RANK == 0)
    {
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
    const int total_steps = end - start;
    const int step_size = total_steps / 50; // Print at 50 percent steps
    int loop_count = 0;

    auto start_t_global_clock = std::chrono::high_resolution_clock::now();

    for(int ii=start; ii<end; ii++)
    {

        auto start_t = std::chrono::high_resolution_clock::now();

        const parameters::FileNames fnames = parameters::get_filenames(ii);

        // Run dynamics START -------------------------------------------------
        if (p.dynamics == "gillespie")
        {
            throw std::runtime_error("gillespie not implemented yet");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (p.dynamics == "standard")
        {
            StandardSimulation standard_sim(fnames, p);
            standard_sim.execute();
        }
        else
        {
            printf("Unsupported dynamics");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Run dynamics END ---------------------------------------------------

        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        std::string dt_string = oss.str();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration_seconds = 
            std::chrono::duration_cast<std::chrono::seconds>(stop - start_t);
        int duration_double_seconds = 
            std::chrono::duration<int>(duration_seconds).count();

        loop_count++;

        if (MPI_RANK == 0)
        {
            if (loop_count % step_size == 0 | loop_count == 1)
            {
                auto stop_g = std::chrono::high_resolution_clock::now();
                auto duration_seconds_g = 
                    std::chrono::duration_cast<std::chrono::seconds>(
                        stop_g - start_t_global_clock);
                int duration_double_seconds_g = 
                    std::chrono::duration<int>(duration_seconds_g).count();
                printf(
                    "%s ~ %s done in %i s (%i/%i) total elapsed %i s\n",
                    dt_string.c_str(), fnames.ii_str.c_str(),
                    duration_double_seconds, loop_count, total_steps, 
                    duration_double_seconds_g
                );
                fflush(stdout);
            }
        }
    }

    MPI_Finalize();
}
