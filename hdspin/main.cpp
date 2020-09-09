#include <fstream>      // std::ofstream
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>

#include "gillespie.h"
#include "standard.h"

int main(int argc, char *argv[])
{

    const std::string target_directory = argv[1];
    const int N_timesteps = atoi(argv[2]);
    const int N_spins = atoi(argv[3]);
    const double beta = atof(argv[4]);
    const double beta_critical = atof(argv[5]);
    const int landscape = atoi(argv[6]);  // 0 for EREM, 1 for REM
    const int dynamics = atoi(argv[7]);  // 0 for standard, 1 for gillespie

    // This will be the minimum job index to resume at
    const int resume_at = atoi(argv[8]);

    // Slurm array task id
    const int slurm_arr_id = atoi(argv[9]);

    // Number of jobs to run on this execution
    const int njobs = atoi(argv[10]);

    // Arguments pertaining to the job itself
    printf("saving data to %s\n", argv[1]);
    printf("N_timesteps = %i\n", N_timesteps);
    printf("N_spins = %i\n", N_spins);
    printf("beta = %.02f\n", beta);
    printf("beta_critical = %.02f\n", beta_critical);
    printf("landscape = %i (0=EREM, 1=REM)\n", landscape);
    printf("dynamics = %i (0=standard, 1=gillespie)\n",
        dynamics);

    std::string ii_str, path, is_path;

    const int start = resume_at + slurm_arr_id * njobs;
    const int end = resume_at + (slurm_arr_id + 1) * njobs;

    printf("job ID's %i -> %i\n", start, end);

    for(int ii = start; ii < end; ii++)
    {
        auto start = std::chrono::high_resolution_clock::now();

        ii_str = std::to_string(ii);
        ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');
        path = target_directory + "/" + ii_str + ".txt";
        is_path = target_directory + "/" + ii_str + "_IS.txt";

        if (dynamics == 1)
        {
            gillespie(path, is_path, N_timesteps, N_spins, beta,
                beta_critical, landscape);
        }
        else if (dynamics == 0)
        {
            standard(path, is_path, N_timesteps, N_spins, beta,
                beta_critical, landscape);
        }
        else
        {
            throw std::runtime_error("Unsupported dynamics flag");
        }

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration_seconds = 
            std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        double duration_double_seconds = 
            std::chrono::duration<double>(duration_seconds).count();

        printf("iter %s done in %.02f s\n", path.c_str(),
            duration_double_seconds);
    }
}
