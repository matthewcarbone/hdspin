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
    const int total_steps = end - start;
    int step_size = total_steps / 50; // Print at 50 percent steps
    if (step_size == 0){step_size = 1;}
    int loop_count = 0;

    auto global_start = std::chrono::high_resolution_clock::now();

    const int starting_seed = p.seed;

    for(int ii=start; ii<end; ii++)
    {

        auto t_start = std::chrono::high_resolution_clock::now();

        const parameters::FileNames fnames = parameters::get_filenames(ii);

        // Change the seed based on the MPI rank, very important for seeded runs!
        // This will be ignored later if p.use_manual_seed is false
        p.seed = starting_seed + ii + MPI_RANK * n_tracers_per_MPI_rank;

        // Run dynamics START -------------------------------------------------
        if (p.dynamics == "gillespie")
        {
            GillespieSimulation gillespie_sim(fnames, p);
            gillespie_sim.execute();
        }
        else if (p.dynamics == "standard")
        {
            StandardSimulation standard_sim(fnames, p);
            standard_sim.execute();
        }
        else if (p.dynamics == "auto")
        {
            // TODO
            printf("Unsupported dynamics (auto)");
            MPI_Abort(MPI_COMM_WORLD, 1);   
        }
        else
        {
            printf("Unsupported dynamics (unknown)");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Run dynamics END ---------------------------------------------------

        const int duration = time_utils::get_time_delta(t_start);

        loop_count++;

        if (MPI_RANK == 0)
        {
            if (loop_count % step_size == 0 | loop_count == 1)
            {
                const std::string dt_string = time_utils::get_datetime();
                const int global_duration = time_utils::get_time_delta(global_start);
                
                printf(
                    "%s ~ %s done in %i s (%i/%i) total elapsed %i s\n", dt_string.c_str(), fnames.ii_str.c_str(), duration, loop_count, total_steps, global_duration
                );
                fflush(stdout);
            }
        }
    }

    MPI_Finalize();
}
