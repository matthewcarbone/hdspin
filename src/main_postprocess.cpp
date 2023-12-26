#include <mpi.h>
#include "processing_utils.h"


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    processing_utils::postprocess();
    MPI_Finalize();
}
