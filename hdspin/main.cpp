#include <fstream>      // std::ofstream
#include <omp.h>
#include "mpi.h"
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
    const int njobs = atoi(argv[8]);
    const int resume_at = atoi(argv[9]);

    int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int iam = 0, np = 1, file_ID;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);

    // Assign the jobs to the appropriate processes
    const int jobs_per_rank = njobs / numprocs;
    const int remainder = njobs % numprocs;
    const int start = rank * jobs_per_rank + resume_at;
    int end = (rank + 1) * jobs_per_rank + resume_at;
    if (rank == numprocs - 1){end += remainder;}

    printf("RANK %i -- program name is %s\n", rank, argv[0]);

    // Arguments pertaining to the job itself
    printf("RANK %i -- saving data to %s\n", rank, argv[1]);
    printf("RANK %i -- N_timesteps = %i\n", rank, N_timesteps);
    printf("RANK %i -- N_spins = %i\n", rank, N_spins);
    printf("RANK %i -- beta = %.02f\n", rank, beta);
    printf("RANK %i -- beta_critical = %.02f\n", rank, beta_critical);
    printf("RANK %i -- landscape = %i (0=EREM, 1=REM)\n", rank, landscape);
    printf("RANK %i -- dynamics = %i (0=standard, 1=gillespie)\n", rank,
        dynamics);

    // The number of jobs total
    printf("RANK %i -- total simulations = %i\n", rank, njobs);

    // We can resume the simulation at a non-zero value
    printf("RANK %i -- resuming at index = %i", rank, resume_at);

    printf("RANK %i -- remainder on this rank = %i\n", rank, remainder);
    printf("RANK %i -- simulations on this rank = %i\n", rank, end - start);
    printf("RANK %i -- start on this rank = %i\n", rank, start);
    printf("RANK %i -- end on this rank = %i\n", rank, end);

    std::string ii_str, path, is_path;

    #pragma omp parallel
    { 
        np = omp_get_num_threads();
        #pragma omp for private(ii_str, path, is_path)
        for(int ii = start; ii < end; ii++)
        {
            iam = omp_get_thread_num();
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

            printf("Thrd %d/%d proc %d/%d on %s iter file %s done in %.02f\n",
                iam, np, rank, numprocs, processor_name, path.c_str(),
                duration_double_seconds);
        }
    }

    MPI_Finalize();
}
