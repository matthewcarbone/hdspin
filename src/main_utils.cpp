#include <chrono>
#include <mpi.h>

#include "utils.h"
#include "spin.h"
#include "obs1.h"
#include "obs2.h"

#include "main_utils.h"


#define DEFAULT_TAG 0
#define MASTER 0
#define DIAGNOSTIC_TAG 1


void step_all_(const double waiting_time, const double simulation_clock, OnePointObservables& obs1, PsiConfig& psiConfig, PsiBasin& psiBasin, AgingConfig& agingConfig, AgingBasin& agingBasin)
{
    obs1.step(waiting_time, simulation_clock);
    psiConfig.step(waiting_time);
    psiBasin.step(waiting_time);
    agingConfig.step(simulation_clock);
    agingBasin.step(simulation_clock);
}


void execute(const int job_index, const utils::SimulationParameters const_params)
{

    const utils::FileNames fnames = utils::get_filenames(job_index);
    utils::SimulationParameters params = const_params;
    const unsigned int starting_seed = const_params.seed;
    params.seed = starting_seed + job_index;

    EnergyMapping emap(params);
    SpinSystem sys(params, emap);

    // Special case of the standard spin dynamics: if rtp.loop_dynamics == 2,
    // then the timestep is divided by rtp.N_spins.
    double waiting_time;

    // Simulation parameters
    double simulation_clock = 0.0;

    OnePointObservables obs1(fnames, params, sys);
    PsiConfig psiConfig(fnames, params, sys);
    PsiBasin psiBasin(fnames, params, sys);
    AgingConfig agingConfig(fnames, params, sys);
    AgingBasin agingBasin(fnames, params, sys);

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

        step_all_(
            waiting_time,
            simulation_clock,
            obs1,
            psiConfig,
            psiBasin,
            agingConfig,
            agingBasin
        );

        if (simulation_clock > params.N_timesteps){break;}
    }
}

double get_sim_time(utils::SimulationParameters p, const std::string dynamics)
{
    double simulation_clock = 0.0;
    p.dynamics = dynamics;
    EnergyMapping emap(p);
    SpinSystem sys(p, emap);
    auto t_start = std::chrono::high_resolution_clock::now();
    for (size_t cc=0; cc<1e7; cc++)
    {
        simulation_clock += sys.step();
        if (simulation_clock > p.N_timesteps){break;}
    }
    return utils::get_time_delta(t_start) / simulation_clock;
}



void master(const size_t min_index_inclusive, const size_t max_index_exclusive)
{

    int mpi_world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);

    const size_t i0 = min_index_inclusive;
    const size_t i1 = max_index_exclusive;
    size_t step_size = (i1 - i0) / 10; // Print at 10 percent steps
    if (step_size == 0){step_size = 1;}

    auto t_start = std::chrono::high_resolution_clock::now();

    int worker_rank, tag;
    size_t loop_count = 0;
    for (size_t next_job=i0; next_job<i1; next_job++)
    {

        loop_count++;

        // Step 2: the master process receives the worker rank and stores it in
        // worker_rank. This indicates to master that worker is ready for a job
        MPI_Recv(&worker_rank, 1, MPI_INT, MPI_ANY_SOURCE, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Step 3: the job is assigned to the worker and the job index is
        // incremented
        MPI_Send(&next_job, 1, MPI_INT, worker_rank, DEFAULT_TAG, MPI_COMM_WORLD);

        if (loop_count % step_size == 0 | loop_count == 1)
        {
            const std::string dt_string = utils::get_datetime();
            const double duration = utils::get_time_delta(t_start);
            std::string ii_str = std::to_string(next_job + 1);
            ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');
            printf(
                "%s ~ %s | %.01f s total\n",
                dt_string.c_str(),
                ii_str.c_str(),
                duration
            );
        }   
    }

    const int END_SIGNAL = -1;

    // We're done on the master process, send termination signals to all of
    // the workers

    std::cout << "!!! DONE !!!" << std::endl;
    std::cout << "Diagnostics:" << std::endl;
    for (size_t cc=0; cc<mpi_world_size - 1; cc++)
    {

        // Receive the ready signal from some source
        MPI_Recv(&worker_rank, 1, MPI_INT, MPI_ANY_SOURCE, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Send the termination signal to that worker, this breaks it out
        // of its core loop
        MPI_Send(&END_SIGNAL, 1, MPI_INT, worker_rank, DEFAULT_TAG, MPI_COMM_WORLD);

        // And finally we ask for the total number of tasks it had to compute
        int tmp_n_jobs;
        MPI_Recv(&tmp_n_jobs, 1, MPI_INT, MPI_ANY_SOURCE, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        std::cout << "RANK " << worker_rank << " finished " << tmp_n_jobs << " jobs" << std::endl;
    }

}


void worker(const utils::SimulationParameters params)
{

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int total_jobs = 0;

    while (true)
    {
        // Step 1: worker `rank` sends its rank to the master process
        MPI_Send(&mpi_rank, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD);

        // Step 4: the job is received from master and executed
        int job_index;
        MPI_Recv(&job_index, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Step 5: if all jobs are completed, the worker receives a termination
        // signal, breaking the process out of the control flow
        if (job_index < 0)
        {
            // Now we can return some information back to the master process about
            // how many jobs were completed.
            MPI_Send(&total_jobs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD);

            // Finish signal received, break from the loop
            return;
        }

        execute(job_index, params);

        total_jobs++;
    }

    
}


// "Public API"

namespace main_utils
{

void auto_determine_dynamics_(utils::SimulationParameters* params)
{
    int mpi_rank, mpi_world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (params->dynamics != "auto"){return;}

    double standard_time = 0.0;
    double gillespie_time = 0.0;
    double standard_std = 0.0;
    double gillespie_std = 0.0;
    const utils::SimulationParameters p = *params;
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
        MPI_Barrier(MPI_COMM_WORLD);

        // Now, we send everything to rank 0
        if (mpi_rank != 0)
        {
            MPI_Send(&times[mpi_rank], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        else
        {
            for (unsigned int rank=1; rank<mpi_world_size; rank++)
            {
                MPI_Recv(&times[rank], 1, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // And let rank 0 deal with the rest
        if (mpi_rank == MASTER)
        {
            std::vector<double> standard_times;
            std::vector<double> gillespie_times;
            for (unsigned int ii=0; ii<mpi_world_size; ii++)
            {
                if (ii % 2 == 0){standard_times.push_back(times[ii]);}
                else{gillespie_times.push_back(times[ii]);}
            }

            // Calculate the mean and standard deviation
            standard_time = utils::mean_vector(standard_times);
            standard_std = sqrt(utils::variance_vector(standard_times));
            gillespie_time = utils::mean_vector(gillespie_times);
            gillespie_std = sqrt(utils::variance_vector(gillespie_times));
        }
    }

    if (mpi_rank == MASTER)
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

    MPI_Bcast(&result_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (result_int == 1){params->dynamics = "gillespie";}
    else if (result_int == 0){params->dynamics = "standard";}
    else{MPI_Abort(MPI_COMM_WORLD, 1);}
}


/**
 * @brief Prints information about the current processor that is running
 * the job on the provided rank
 * 
 * @param mpi_rank
 * @param mpi_world_size
 */
// void print_processor_information(const int mpi_rank, const int mpi_world_size)
// {
//     // Get the name of the processor
//     char processor_name[MPI_MAX_PROCESSOR_NAME];
//     int name_len;

//     // Replace name_len with NULL?
//     MPI_Get_processor_name(processor_name, &name_len);

//     // Print off a hello world message
//     printf("Ready: processor %s, rank %d/%d\n", processor_name, mpi_rank, mpi_world_size);
// }


void initialize_grids_and_config(const utils::SimulationParameters p)
{

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (mpi_rank == MASTER)
    {
        utils::log_parameters(p);
        utils::make_directories();
        utils::make_energy_grid_logspace(p.log10_N_timesteps, p.grid_size);
        utils::make_pi_grids(p.log10_N_timesteps, p.dw, p.grid_size);
        json j = utils::parameters_to_json(p);
        std::ofstream o("config.json");
        o << std::setw(4) << j << std::endl;
    }

    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
}

void execute_process_pool(const utils::SimulationParameters params)
{
    int mpi_rank, mpi_world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);

    // Execute the process pool
    if (mpi_rank == MASTER)
    {
        // TODO add some logic for checkpoint-restart here
        std::cout << "Running " << mpi_world_size << " procs. " << mpi_world_size - 1 << " compute, 1 controller " << std::endl;
        master(0, params.n_tracers);
    }
    else
    {
        worker(params);
    }
}

}
