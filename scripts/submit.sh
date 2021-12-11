#!/bin/bash

mkdir job_data
jobID=$(sbatch calc.sbatch.sh)
sbatch --dependency=afterany:"${jobID##* }" eval_1.sbatch.sh
