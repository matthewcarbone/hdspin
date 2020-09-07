#include <fstream>      // std::ofstream
#include <omp.h>
#include "mpi.h"
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>

#include "gillespie.h"

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

    printf("program name is %s\n", argv[0]);

    // Arguments pertaining to the job itself
    printf("saving data to %s\n", argv[1]);
    printf("N_timesteps = %i\n", N_timesteps);
    printf("N_spins = %i\n", N_spins);
    printf("beta = %.02f\n", beta);
    printf("beta_critical = %.02f\n", beta_critical);
    printf("landscape = %i (0=EREM, 1=REM)\n", landscape);
    printf("dynamics = %i (0=standard, 1=gillespie)\n", dynamics);

    // The number of jobs total
    printf("total simulations = %i\n", njobs);

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
    const int start = rank * jobs_per_rank;
    int end = (rank + 1) * jobs_per_rank;
    if (rank == numprocs - 1){end += remainder;}

    printf("remainder on this rank = %i\n", remainder);
    printf("simulations on this rank = %i\n", end - start);
    printf("start on this rank = %i\n", start);
    printf("end on this rank = %i\n", end);

    std::string ii_str, path;

    #pragma omp parallel
    { 
        np = omp_get_num_threads();
        #pragma omp for private(ii_str, path)
        for(int ii = start; ii < end; ii++)
        {
            iam = omp_get_thread_num();
            auto start = std::chrono::high_resolution_clock::now();

            ii_str = std::to_string(ii);
            ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');
            path = target_directory + "/" + ii_str + ".txt";

            if (dynamics == 1)
            {
                gillespie(path, N_timesteps, N_spins, beta, beta_critical,
                    landscape);
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
