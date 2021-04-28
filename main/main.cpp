#include <fstream>      // std::ofstream
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <sstream>
#include <assert.h>
#include <mpi.h>

#include "Utils/structures.h"
#include "Engine/sim.h"


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
    MPI_Barrier(MPI_COMM_WORLD);

    // Parse all the input
    const std::string target_directory = argv[1];
    const std::string grids_directory = argv[2];
    const std::string final_directory = argv[3];
    const int log_N_timesteps = atoi(argv[4]);
    const int N_spins = atoi(argv[5]);
    const double beta = atof(argv[6]);
    const double beta_critical = atof(argv[7]);
    const std::string str_landscape = argv[8];  // 0 for EREM, 1 for REM
    const std::string str_dynamics = argv[9];  // 0 for standard, 1 for gillespie
    const int n_tracers = atoi(argv[10]);

    int landscape = -1;
    if (str_landscape == "erem"){landscape = 0;}
    else if (str_landscape == "rem"){landscape = 1;}
    assert(landscape > -1);

    // Indexes the Standard (0) or Gillespie (1) protocols.
    int dynamics = -1;

    // 0 for standard, 1 for loop over N, 2 for standard dynamics but where the
    // timestep is divided by N.
    int loop_dynamics = -1;
    if (str_dynamics == "standard")
    {
        dynamics = 0;
        loop_dynamics = 0;
    }
    else if (str_dynamics == "standard-loop")
    {
        dynamics = 0;
        loop_dynamics = 1;
    }
    else if (str_dynamics == "standard-divN")
    {
        dynamics = 0;
        loop_dynamics = 2;
    }
    else if (str_dynamics == "gillespie")
    {
        dynamics = 1;
        loop_dynamics = 0;
    }
    else if (str_dynamics == "gillespie-divN")
    {
        dynamics = 1;
        loop_dynamics = 2;
    }
    assert(dynamics > -1);
    assert(loop_dynamics > -1);

    const RuntimeParameters params = get_runtime_params(log_N_timesteps,
        N_spins, beta, beta_critical, landscape, loop_dynamics);

    // Arguments pertaining to the job itself
    if (MPI_RANK == 0)
    {
        printf("* Saving data to: %s\n", argv[1]);
        printf("* log10(N) timesteps: %i\n", params.log_N_timesteps);
        printf("* N timesteps: %lli\n", params.N_timesteps);
        printf("* N spins: %i\n", params.N_spins);
        printf("* N configs: %lli\n", params.N_configs);
        printf("* 1/T: c.a. %.02f\n", params.beta);
        printf("* 1/Tc: c.a. %.02f\n", params.beta_critical);
        printf("* Landscape: %s (%i)\n", str_landscape.c_str(),
            params.landscape);
        printf("* Dynamics: %s (%i)\n", str_dynamics.c_str(), dynamics);
        printf("* Loop dynamics: %i\n", loop_dynamics);
        printf("* Energetic threshold: c.a. %.05e\n",
            params.energetic_threshold);
        printf("* Entropic attractor: c.a. %.05e\n",
            params.entropic_attractor);

        if ((landscape == 0) &&
            (beta >= 2.0 * beta_critical | beta <= beta_critical))
        {
            printf(
                "WARNING: beta restriction bc < b < 2bc not satisfied,\n"
            );
            printf("WARNING: => entropic attractor set arbitrarily high\n");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);  // So everything prints cleanly

    const int resume_at = 0;
    const int start = resume_at + MPI_RANK * n_tracers;
    const int end = resume_at + (MPI_RANK + 1) * n_tracers;

    printf("RANK %i/%i job ID's %i -> %i\n", MPI_RANK, MPI_WORLD_SIZE, start,
        end);

    MPI_Barrier(MPI_COMM_WORLD);

    // Define some helpers to be used to track progress.
    const int total_steps = end - start;
    const int step_size = total_steps / 50; // Print at 50 percent steps
    int loop_count = 0;

    auto start_t_global_clock = std::chrono::high_resolution_clock::now();

    for(int ii = start; ii < end; ii++)
    {
        auto start_t = std::chrono::high_resolution_clock::now();

        const FileNames fnames = get_filenames(ii, target_directory,
            grids_directory);

        // Run dynamics START -------------------------------------------------
        if (dynamics == 1)
        {
            GillespieSimulation gillespie_sim(fnames, params);
            gillespie_sim.execute();
        }
        else if (dynamics == 0)
        {
            StandardSimulation standard_sim(fnames, params);
            standard_sim.execute();
        }
        else
        {
            printf("Unsupported dynamics flag");
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
