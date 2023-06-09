# hdspin
REM/EREM sandbox

## Installation instructions

The `hdspin` repository requires no external libraries whatsoever. The file `json.hpp` is now included in the repository as per the terms of the MIT license at `inc/Json/json.hpp`.

To make on Linux should be as simple as

```bash
make local_linux
```


## How to run

### Local
On your local machine, running `hdspin` is very straightforward. Create a directory, e.g., `scratch/my_test_hdspin_run`, and enter it.

Inside this directory there need be only a single file, `input.json`. This file must contain _all_ of the following keywords in `json` format (basically a python dictionary). For example,

```json
{
    "log_N_timesteps": 7,
    "N_spins": 20,
    "beta": 1.667,
    "landscape" : "REM-num",
    "dynamics": "gillespie",
    "divN": false,
    "memory": -1,
    "max_ridges": 200,
    "n_tracers": 100,
    "memoryless_retain_last_energy": false
}
```

The above tells the code to do the following (this should clarify what each key is):

* Run for `10^7` time steps (can be any positive number)
* Use a spin system of `20` spins
* Use inverse temperature `beta/beta_c = 1.667`. The provided value of beta is always actually in units of the critical beta.
    * `beta_c = 1` for EREM
    * `beta_c = sqrt(2 ln 2) ~ 1.177`
* Use the Random Energy Model landscape (Gaussian energy landscape) but with calculating the REM threshold energy numerically with 10k (not an input parameter, this is hard-coded) samples to compute it on average. Options are `REM`, `REM-num` and `EREM`.
* Use Gillespie dynamics. Options are `standard` and `gillespie`
* Do not use divN dynamics. Timesteps are standard. If `true`, will divide the waiting time by N.
* Use a full memory simulation.
    * If `memory` is `-1`, this means to initialize the energy landscape at instantiation, which remains fixed for the simulation duration. This will require `O(2^N)` memory however, so it cannot be used when `N` is large.
    * If `memory` is `0`, this means to run a memoryless simulation. Inherent structure observables will not be calculated, and the energy is not state-dependent.
    * If `memory` is greater than 0, this is a special case where we use a least-recently used queue to save only `memory` samples in the energy landscape. For example, if `memory` is 100, then only the 100 most recently visited configurations' energies will be kept. If a new configuration is visited, one that is not in the queue, it will be randomly resampled, even if it was visited a long time ago.
* Compute `200` max ridges only, then stop saving.
* Use `100` tracers _per MPI rank_!!
* Do not use the special case of the memory dynamics where we retain only the last-visited neighbor in the neighbor list.

It is often convenient to have a submit script as well. For local jobs, all we need is

```bash
#!/bin/bash

exe_path=/Users/mc/Github/hdspin/exe  # <- For example

# Execution phase
mpiexec -n 4 "$exe_path"/main.out > output.log 2> warn.log

# Evaluation phase
python3 $exe_path/eval.py >> output.log 2>> warn.log

if [[ -e "data/00000000_ridge_E.txt" ]]; then
    cat data/*_ridge_E.txt >> final/ridge_E.txt
fi

if [[ -e "data/00000000_ridge_S.txt" ]]; then
    cat data/*_ridge_S.txt >> final/ridge_S.txt
fi
````

This launches 4 MPI processes (which are basically embarrassingly parallel), each of which will have `n_tracers` tracers as specified in the input file. Once that is complete, there is another script, `eval.py` which will run postprocessing on most of the data.


### HPC

On a high-performance computing platform, the submit script is just a bit more complicated, including e.g. the SLURM runtime magic. One should run the `submit.sh` file, which is shown just below. Below that, are the two other helper scripts which are submitted as dependencies.

**`submit.sh`**
```bash
#!/bin/bash

mkdir job_data
jobID=$(sbatch calc.sbatch.sh)
sbatch --dependency=afterany:"${jobID##* }" eval_1.sbatch.sh
```

**`calc.sbatch.sh`**
```bash
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
#SBATCH --exclude=compute-0-6

exe_path=/home/mcarbone/hdspin/exe

module load compilers/gcc-10.1.0
export OMP_NUM_THREADS=1

mpiexec "$exe_path"/main.out > output.log 2> warn.log
```

**`eval_1.sbatch.sh`**
```bash
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
source activate py37mpi  # Don't forget to load packages if needed

python3 $exe_path/eval.py >> output.log 2>> warn.log

if [[ -e "data/00000000_ridge_E.txt" ]]; then
    cat data/*_ridge_E.txt >> final/ridge_E.txt
fi

if [[ -e "data/00000000_ridge_S.txt" ]]; then
    cat data/*_ridge_S.txt >> final/ridge_S.txt
fi
```





























