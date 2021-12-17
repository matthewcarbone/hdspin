#!/bin/bash

#SBATCH --nice=10000
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -n 1
#SBATCH -p New
#SBATCH --mem=7GB
#SBATCH --job-name=eval
#SBATCH --output=job_data/eval_%A.out
#SBATCH --error=job_data/eval_%A.err
#SBATCH --exclude=compute-0-6

exe_path=/home/mcarbone/hdspin/exe
source activate py37mpi

python3 $exe_path/eval.py >> output.log 2>> warn.log

if [[ -e "data/00000000_ridge_E.txt" ]]; then
    cat data/*_ridge_E.txt >> final/ridge_E.txt
fi

if [[ -e "data/00000000_ridge_S.txt" ]]; then
    cat data/*_ridge_S.txt >> final/ridge_S.txt
fi
