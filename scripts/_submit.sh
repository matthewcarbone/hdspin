#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p New
#SBATCH --mem=2GB
#SBATCH --output=job_data/rem.%A.out
#SBATCH --error=job_data/rem.%A.err
#SBATCH --job-name=rem

# export OMP_NUM_THREADS=1
# export USE_SIMPLE_THREADED_LEVEL3=1

#export PATH="/opt/pb/gcc-10.1.0/bin:$PATH"
#export LD_LIBRARY_PATH="/opt/pb/gcc-10.1.0/lib:/opt/pb/gcc-10.1.0/lib64:$LD_LIBRARY_PATH"

module load compilers/gcc-10.1.0
./main.o "$@"
