#!/bin/bash

#SBATCH --nice=10000
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -n 24
#SBATCH -p New
#SBATCH --mem=62GB
#SBATCH --job-name=calc
#SBATCH --output=job_data/calc_%A.out
#SBATCH --error=job_data/calc_%A.err

exe_path=/home/mcarbone/hdspin/exe

module load compilers/gcc-10.1.0
export OMP_NUM_THREADS=1

mpiexec "$exe_path"/main.out > output.log 2> warn.log
