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

// Requires nlohmann::json
// See here: https://github.com/nlohmann/json
// Download the json.hpp file and put it in your path, it's that easy
#include <json.hpp>
using json = nlohmann::json;


RuntimeParameters get_runtime_parameters()
{

    std::ifstream ifs("input.json");
    json inp = json::parse(ifs);

    RuntimeParameters rtp;

    // GEneral
    rtp.log_N_timesteps = int(inp["log_N_timesteps"]);
    rtp.N_timesteps = ipow(10, rtp.log_N_timesteps);
    rtp.N_spins = int(inp["N_spins"]);
    rtp.beta = float(inp["beta"]);

    // Landscape
    rtp.landscape = inp["landscape"];
    if (rtp.landscape != "EREM" && rtp.landscape != "EREM")
    {
        throw std::runtime_error("Invalid landscape");
    }

    // Beta critical
    if (rtp.landscape == 0){rtp.beta_critical = 1.0;}
    else{rtp.beta_critical = 1.177410022515475;}  // This is ~sqrt(2 ln 2)

    // Dynamics
    std::string _dynamics = inp["dynamics"];
    if (_dynamics == "standard"){rtp.dynamics_flag = 0;}
    else if (_dynamics == "gillespie"){rtp.dynamics_flag = 1;}
    else{throw std::runtime_error("Invalid dynamics");}

    // Divide-by-N status
    rtp.divN = int(inp["divN"]);
    if (rtp.divN > 1 || rtp.divN < 0){throw std::runtime_error("Invalid divN");}

    // Memory status
    rtp.memory = int(inp["memory"]);
    if (rtp.memory < -1){throw std::runtime_error("Invalid memory");}

    if (rtp.memory != 0){rtp.N_configs = ipow(2, rtp.N_spins);}
    else{rtp.N_configs = -1;}  // Don't use in the memoryless case

    // Maximum number of ridges before recording stops
    rtp.max_ridges = int(inp["max_ridges"]);
    if (rtp.max_ridges < -1){throw std::runtime_error("Invalid max_ridges");}

    // Get the energy barrier information
    double et, ea;
    rtp.valid_entropic_attractor = true;
    if (rtp.landscape == "EREM") // EREM
    {
        et = -1.0 / rtp.beta_critical * log(rtp.N_spins);
        ea = 1.0 / (rtp.beta - rtp.beta_critical)
            * log((2.0 * rtp.beta_critical - rtp.beta) / rtp.beta_critical);

        if (rtp.beta >= 2.0 * rtp.beta_critical | rtp.beta <= rtp.beta_critical)
        {
            ea = 1e16;  // Set purposefully invalid value instead of nan or inf
            rtp.valid_entropic_attractor = false;
        }
    }
    else // REM
    {
        et = -sqrt(2.0 * rtp.N_spins * log(rtp.N_spins));
        ea = -rtp.N_spins * rtp.beta / 2.0;
    }

    rtp.energetic_threshold = et;
    rtp.entropic_attractor = ea;

    rtp.n_tracers = inp["n_tracers"];

    return rtp;
}

void log_rtp(const RuntimeParameters rtp)
{
    printf("log_N_timesteps = %i\n", rtp.log_N_timesteps);
    printf("N_timesteps = %lli\n", rtp.N_timesteps);
    printf("N_spins = %i\n", rtp.N_spins);
    printf("N_configs = %lli\n", rtp.N_configs);
    printf("beta = %.05f\n", rtp.beta);
    printf("beta_critical = %.05f\n", rtp.beta_critical);

    printf("landscape = %s\n", rtp.landscape.c_str());
    std::string _dynamics;

    if (rtp.dynamics_flag == 0){_dynamics = "standard";}
    else if (rtp.dynamics_flag == 1){_dynamics = "gillespie";}
    printf("dynamics = %i (%s)\n", rtp.dynamics_flag, _dynamics.c_str());

    printf("divN = %i\n", rtp.divN);

    std::string _mem;
    if (rtp.memory == 0){_mem = "memoryless simulation";}
    else if (rtp.memory == -1){_mem = "full initialization of memory";}
    else{_mem = "specific memory value set";}
    printf("memory = %i (%s)\n", rtp.memory, _mem.c_str());
    printf("max E ridges = %i\n", rtp.max_ridges);
    printf("energetic threshold = %.05f\n", rtp.energetic_threshold);
    printf("entropic attractor = %.05f\n", rtp.entropic_attractor);
    printf("valid_entropic_attractor = %i\n", rtp.valid_entropic_attractor);
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

    RuntimeParameters rtp = get_runtime_parameters();

    const int n_tracers = rtp.n_tracers;

    // Get the information for this MPI rank
    const int resume_at = 0;
    const int start = resume_at + MPI_RANK * n_tracers;
    const int end = resume_at + (MPI_RANK + 1) * n_tracers;

    // So everything prints cleanly
    printf("RANK %i/%i ID's %i -> %i\n", MPI_RANK, MPI_WORLD_SIZE, start, end);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_RANK == 0){log_rtp(rtp);}
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    exit(0);

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
        if (rtp.dynamics_flag % 2 == 1)
        {
            GillespieSimulation gillespie_sim(fnames, rtp);
            gillespie_sim.execute();
        }
        else if (rtp.dynamics_flag % 2 == 0)
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
