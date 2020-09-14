#include <fstream>      // std::ofstream
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <sstream>

#include "gillespie.h"
#include "standard.h"
#include "utils/grid_utils.h"


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

    // This will be the minimum job index to resume at
    const int resume_at = atoi(argv[9]);

    // Slurm array task id
    const int slurm_arr_id = atoi(argv[10]);

    // Number of jobs to run on this execution
    const int njobs = atoi(argv[11]);

    // Arguments pertaining to the job itself
    printf("saving data to %s\n", argv[1]);
    printf("log_N_timesteps = %i\n", log_N_timesteps);
    printf("N_spins = %i\n", N_spins);
    printf("beta = %.02f\n", beta);
    printf("beta_critical = %.02f\n", beta_critical);
    printf("landscape = %i (0=EREM, 1=REM)\n", landscape);
    printf("dynamics = %i (0=standard, 1=gillespie)\n", dynamics);

    std::string ii_str, e_path, psi_config_path;

    const int start = resume_at + slurm_arr_id * njobs;
    const int end = resume_at + (slurm_arr_id + 1) * njobs;

    printf("job ID's %i -> %i\n", start, end);
    for(int ii = start; ii < end; ii++)
    {
        auto start = std::chrono::high_resolution_clock::now();

        ii_str = std::to_string(ii);
        ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');

        // Define all the observable trackers ---------------------------------
        // Energy
        e_path = target_directory + "/" + ii_str + "_energy.txt";
        EnergyGrid energy_grid(grids_directory + "/energy.txt");
        energy_grid.open_outfile(e_path);
        // Psi config
        psi_config_path = target_directory + "/" + ii_str + "_psi_config.txt";
        PsiConfigCounter psi_config_counter(log_N_timesteps, psi_config_path);

        // --------------------------------------------------------------------

        if (dynamics == 1)
        {
            // printf("Running Gillespie dynamics\n");
            gillespie(energy_grid, psi_config_counter, log_N_timesteps,
                N_spins, beta, beta_critical, landscape);
        }

        else if (dynamics == 0)
        {
            // printf("Running standard dynamics\n");
            standard(energy_grid, psi_config_counter, log_N_timesteps, N_spins,
                beta, beta_critical, landscape);
        }

        else
        {
            throw std::runtime_error("Unsupported dynamics flag");
        }

        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        std::string dt_string = oss.str();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration_seconds = 
            std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        double duration_double_seconds = 
            std::chrono::duration<double>(duration_seconds).count();

        // Close the outfiles and write to disk when not doing so dynamically
        energy_grid.close_outfile();
        psi_config_counter.write_to_disk();

        printf("%s ~ %s done in %.02f s\n", dt_string.c_str(), ii_str.c_str(),
            duration_double_seconds);
    }
}
