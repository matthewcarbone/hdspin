#include <fstream>      // std::ofstream
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <sstream>
#include <omp.h>

#include "Utils/structures.h"
#include "Engine/sim.h"


int main(int argc, char *argv[])
{
    const std::string target_directory = argv[1];
    const std::string grids_directory = argv[2];
    const int log_N_timesteps = atoi(argv[3]);
    const int N_spins = atoi(argv[4]);
    const double beta = atof(argv[5]);
    const double beta_critical = atof(argv[6]);
    const int landscape = atoi(argv[7]);  // 0 for EREM, 1 for REM
    const int dynamics = atoi(argv[8]);  // 0 for standard, 1 for gillespie

    // 0 for standard, 1 for loop over N
    const int loop_dynamics = atoi(argv[9]);
    const RuntimeParameters params = get_runtime_params(log_N_timesteps,
        N_spins, beta, beta_critical, landscape, loop_dynamics);

    // This will be the minimum job index to resume at
    const int resume_at = atoi(argv[10]);

    // Slurm array task id
    const int slurm_arr_id = atoi(argv[11]);

    // Number of jobs to run on this execution
    const int njobs = atoi(argv[12]);

    // Arguments pertaining to the job itself
    printf("saving data to %s\n", argv[1]);
    printf("log_N_timesteps = %i\n", log_N_timesteps);
    printf("N_spins = %i\n", N_spins);
    printf("beta = %.02f\n", beta);
    printf("beta_critical = %.02f\n", beta_critical);
    printf("landscape = %i (0=EREM, 1=REM)\n", landscape);
    printf("dynamics = %i (0=standard, 1=gillespie)\n", dynamics);
    printf("loopN = %i\n", loop_dynamics);

    const int start = resume_at + slurm_arr_id * njobs;
    const int end = resume_at + (slurm_arr_id + 1) * njobs;

    printf("job ID's %i -> %i\n", start, end);

    // Define some helpers to be used to track progress.
    const int total_steps = end - start;
    const int step_size = total_steps / 50; // Print at 50 percent steps
    int loop_count = 0;

    auto start_t_global_clock = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
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
        else{throw std::runtime_error("Unsupported dynamics flag");}
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

        #pragma omp atomic
        loop_count++;

        if (loop_count % step_size == 0 | loop_count == 1)
        {
            #pragma omp critical
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
}
