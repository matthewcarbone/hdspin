#include <fstream>      // std::ofstream
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <sstream>
#include <assert.h>
#include <mpi.h>

#include "Engine/sim.h"
#include "Utils/structures.h"
#include "Utils/utils.h"


RuntimeParameters get_runtime_parameters(char *argv[])
{
    RuntimeParameters rtp;

    // Parse all the inputs
    rtp.log_N_timesteps = atoi(argv[1]);
    rtp.N_timesteps = ipow(10, rtp.log_N_timesteps);
    rtp.N_spins = atoi(argv[2]);
    rtp.N_configs = ipow(2, rtp.N_spins);
    rtp.beta = atof(argv[3]);
    rtp.beta_critical = atof(argv[4]);

    // 0 for EREM, 1 for REM
    rtp.landscape = atoi(argv[5]);
    assert(rtp.landscape > -1);
    assert(rtp.landscape < 2);

    // Dynamics flag:
    // standard = 0
    // gillespie = 1
    // standard divN = 2
    // gillespie divN = 3
    rtp.dynamics_flag = atoi(argv[6]);
    assert(rtp.dynamics_flag > -1);
    assert(rtp.dynamics_flag < 4);

    // standard = 0
    // memoryless = 1
    rtp.memoryless = atoi(argv[7]);
    rtp.max_ridges = atoi(argv[8]);

    // Get the energy barrier information
    double et, ea;
    if (rtp.landscape == 0) // EREM
    {
        et = -1.0 / rtp.beta_critical * log(rtp.N_spins);
        ea = 1.0 / (rtp.beta - rtp.beta_critical)
            * log((2.0 * rtp.beta_critical - rtp.beta) / rtp.beta_critical);

        if (rtp.beta >= 2.0 * rtp.beta_critical | rtp.beta <= rtp.beta_critical)
        {
            ea = 1e16;  // Set purposefully invalid value instead of nan or inf
        }
    }
    else if (rtp.landscape == 1) // REM
    {
        et = -sqrt(2.0 * rtp.N_spins * log(rtp.N_spins));
        ea = -rtp.N_spins * rtp.beta / 2.0;
    }
    else
    {
        throw std::runtime_error("Unknown landscape flag");
    }

    rtp.energetic_threshold = et;
    rtp.entropic_attractor = ea;

    return rtp;
}


FileNames get_filenames(const int ii)
{
    std::string ii_str = std::to_string(ii);
    ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');
    FileNames fnames;
    fnames.energy = "data/" + ii_str + "_energy.txt";
    fnames.psi_config = "data/" + ii_str + "_psi_config.txt";
    fnames.psi_config_IS = "data/" + ii_str + "_psi_config_IS.txt";

    fnames.psi_basin_E = "data/" + ii_str + "_psi_basin_E.txt";
    fnames.psi_basin_S = "data/" + ii_str + "_psi_basin_S.txt";
    fnames.psi_basin_E_IS = "data/" + ii_str + "_psi_basin_E_IS.txt";
    fnames.psi_basin_S_IS = "data/" + ii_str + "_psi_basin_S_IS.txt";

    fnames.aging_config_1 = "data/" + ii_str + "_pi1_config.txt";
    fnames.aging_config_2 = "data/" + ii_str + "_pi2_config.txt";

    fnames.aging_basin_1 = "data/" + ii_str + "_pi1_basin.txt";
    fnames.aging_basin_2 = "data/" + ii_str + "_pi2_basin.txt";

    fnames.ridge_E = "data/" + ii_str + "_ridge_E.txt";
    fnames.ridge_S = "data/" + ii_str + "_ridge_S.txt";
    fnames.ridge_E_all = "data/" + ii_str + "_ridge_E_all.txt";
    fnames.ridge_S_all = "data/" + ii_str + "_ridge_S_all.txt";

    fnames.ridge_E_IS = "data/" + ii_str + "_ridge_E_IS.txt";
    fnames.ridge_S_IS = "data/" + ii_str + "_ridge_S_IS.txt";
    fnames.ridge_E_proxy_IS = "data/" + ii_str + "_ridge_E_proxy_IS.txt";
    fnames.ridge_S_proxy_IS = "data/" + ii_str + "_ridge_S_proxy_IS.txt";

    fnames.ii_str = ii_str;
    fnames.grids_directory = "grids";
    return fnames;
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
    printf("Ready: processor %s, rank %d/%d\n", processor_name, MPI_RANK,
        MPI_WORLD_SIZE);

    // Quick barrier to make sure the printing works out cleanly
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    RuntimeParameters rtp = get_runtime_parameters(argv);
    const int n_tracers = atoi(argv[9]);

    // Get the information for this MPI rank
    const int resume_at = 0;
    const int start = resume_at + MPI_RANK * n_tracers;
    const int end = resume_at + (MPI_RANK + 1) * n_tracers;

    printf("RANK %i/%i ID's %i -> %i\n", MPI_RANK, MPI_WORLD_SIZE, start, end);

    // So everything prints cleanly
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // Define some helpers to be used to track progress.
    const int total_steps = end - start;
    const int step_size = total_steps / 50; // Print at 50 percent steps
    int loop_count = 0;

    auto start_t_global_clock = std::chrono::high_resolution_clock::now();

    for(int ii = start; ii < end; ii++)
    {
        auto start_t = std::chrono::high_resolution_clock::now();

        const FileNames fnames = get_filenames(ii);

        // Run dynamics START -------------------------------------------------
        if (rtp.dynamics_flag  == 1)
        {
            GillespieSimulation gillespie_sim(fnames, rtp);
            gillespie_sim.execute();
        }
        else if (rtp.dynamics_flag  == 0)
        {
            StandardSimulation standard_sim(fnames, rtp);
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
