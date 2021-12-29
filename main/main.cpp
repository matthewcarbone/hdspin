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

    // General
    rtp.log_N_timesteps = int(inp["log_N_timesteps"]);
    rtp.N_timesteps = ipow(10, rtp.log_N_timesteps);
    rtp.N_spins = int(inp["N_spins"]);

    // Landscape
    rtp.landscape = inp["landscape"];
    if (rtp.landscape != "EREM" && rtp.landscape != "REM" && rtp.landscape != "REM-num")
    {
        throw std::runtime_error("Invalid landscape");
    }

    // Beta critical
    if (rtp.landscape == "EREM"){rtp.beta_critical = 1.0;}
    else{rtp.beta_critical = 1.177410022515475;}  // This is ~sqrt(2 ln 2)

    // The provided beta is actually beta/betac
    rtp.beta = float(inp["beta"]) * rtp.beta_critical;

    // Dynamics
    rtp.dynamics = inp["dynamics"];
    if (rtp.dynamics != "standard" && rtp.dynamics != "gillespie")
    {
        throw std::runtime_error("Invalid dynamics");
    }

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

        // Find the threshold energy for the REM model numerically. This is
        // more accurate at N < inf than using the analytic version
        if (rtp.landscape == "REM-num")
        {
            std::mt19937 gen;
            std::normal_distribution<double> dist;

            const double p = sqrt(rtp.N_spins);
            dist.param(std::normal_distribution<double>::param_type(0.0, p));
            const double nreps = 10000;
            std::vector<double> vec(rtp.N_spins, 0.0);
            double min_val;
            double tmp = 0.0;
            for (int rep=0; rep<nreps; rep++)
            {
                for (int nn=0; nn<rtp.N_spins; nn++)
                {
                    vec[nn] = dist(gen);
                }
                min_val = *std::min_element(vec.begin(), vec.end());
                tmp += min_val;
            }
            et = tmp / ((double) nreps);
        }

        // Else just find it analytically
        else{et = -sqrt(2.0 * rtp.N_spins * log(rtp.N_spins));}
        ea = -rtp.N_spins * rtp.beta / 2.0;
    }

    rtp.energetic_threshold = et;
    rtp.entropic_attractor = ea;

    rtp.n_tracers = int(inp["n_tracers"]);

    rtp.memoryless_retain_last_energy = int(inp["memoryless_retain_last_energy"]);

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
    printf("dynamics = %s\n", rtp.dynamics.c_str());

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
    printf("if memoryless, retaining last energy = %i\n", rtp.memoryless_retain_last_energy);
}


FileNames get_filenames(const int ii)
{
    std::string ii_str = std::to_string(ii);
    ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');

    FileNames fnames;

    fnames.energy = "data/" + ii_str + "_energy.txt";
    fnames.energy_avg_neighbors = "data/" + ii_str + "_energy_avg_neighbors.txt";
    fnames.energy_IS = "data/" + ii_str + "_energy_IS.txt";

    fnames.psi_config = "data/" + ii_str + "_psi_config.txt";
    fnames.psi_config_IS = "data/" + ii_str + "_psi_config_IS.txt";

    fnames.psi_basin_E = "data/" + ii_str + "_psi_basin_E.txt";
    fnames.psi_basin_S = "data/" + ii_str + "_psi_basin_S.txt";
    fnames.psi_basin_E_IS = "data/" + ii_str + "_psi_basin_E_IS.txt";
    fnames.psi_basin_S_IS = "data/" + ii_str + "_psi_basin_S_IS.txt";

    fnames.aging_config = "data/" + ii_str + "_aging_config.txt";
    fnames.aging_config_IS = "data/" + ii_str + "_aging_config_IS.txt";

    fnames.aging_basin_E = "data/" + ii_str + "_aging_basin_E.txt";
    fnames.aging_basin_E_IS = "data/" + ii_str + "_aging_basin_E_IS.txt";
    fnames.aging_basin_S = "data/" + ii_str + "_aging_basin_S.txt";
    fnames.aging_basin_S_IS = "data/" + ii_str + "_aging_basin_S_IS.txt";

    fnames.ridge_E_all = "data/" + ii_str + "_ridge_E.txt";
    fnames.ridge_S_all = "data/" + ii_str + "_ridge_S.txt";
    fnames.ridge_E_IS_all = "data/" + ii_str + "_ridge_E_IS.txt";
    fnames.ridge_S_IS_all = "data/" + ii_str + "_ridge_S_IS.txt";

    fnames.unique_configs = "data/" + ii_str + "_unique_configs.txt";

    fnames.ii_str = ii_str;
    fnames.grids_directory = "grids";
    return fnames;
}


void make_directories()
{
    system("mkdir data");
    system("mkdir grids");
}

void make_energy_grid_logspace(const int log10_timesteps,
    const int n_gridpoints)
{
    std::vector <long long> v;
    v.push_back(0.0);
    const double delta = ((double) log10_timesteps) / ((double) n_gridpoints);
    for (int ii=0; ii<n_gridpoints + 1; ii++)
    {
        int val = int(pow(10, ((double) ii * delta)));
        v.push_back(val);
    }

    v.erase(unique(v.begin(), v.end()), v.end());

    const std::string fname = "grids/energy.txt";
    FILE* outfile = fopen(fname.c_str(), "w");
    for (int ii=0; ii<v.size(); ii++)
    {
        fprintf(outfile, "%lli\n", v[ii]);
    }
    fclose(outfile);
}

void make_pi_grids(const int log10_timesteps, const double dw,
    const int n_gridpoints)
{
    std::vector <long long> v1;
    std::vector <long long> v2;
    const int nMC = int(pow(10, log10_timesteps));
    const int tw_max = int(nMC / (dw + 1.0));
    const double delta = ((double) log10(tw_max)) / ((double) n_gridpoints);

    int _v1;
    for (int ii=0; ii<n_gridpoints + 1; ii++)
    {
        _v1 = int(pow(10, ((double) ii * delta)));
        v1.push_back(_v1);
    }
    v1.erase(unique(v1.begin(), v1.end()), v1.end());

    int _v2;
    for (int ii=0; ii<v1.size(); ii++)
    {
        _v2 = int(v1[ii] * (dw + 1.0));
        v2.push_back(_v2);
    }

    const std::string fname1 = "grids/pi1.txt";
    const std::string fname2 = "grids/pi2.txt";
    FILE* outfile1 = fopen(fname1.c_str(), "w");
    FILE* outfile2 = fopen(fname2.c_str(), "w");

    for (int ii=0; ii<v1.size(); ii++)
    {
        fprintf(outfile1, "%lli\n", v1[ii]);
        fprintf(outfile2, "%lli\n", v2[ii]);
    }

    fclose(outfile1);
    fclose(outfile2);
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

    if (MPI_RANK == 0)
    {
        printf("----------------------------------------------------------\n");
        log_rtp(rtp);
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_RANK == 0)
    {
        make_directories();
        make_energy_grid_logspace(rtp.log_N_timesteps, 100);
        make_pi_grids(rtp.log_N_timesteps, 0.5, 100);
        printf("----------------------------------------------------------\n");
    }
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
        if (rtp.dynamics == "gillespie")
        {
            GillespieSimulation gillespie_sim(fnames, rtp);
            gillespie_sim.execute();
        }
        else if (rtp.dynamics == "standard")
        {
            StandardSimulation standard_sim(fnames, rtp);
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
