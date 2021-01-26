#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import numpy as np

import os
import subprocess
import yaml
import warnings


MAX_LONG = 9223372036854775807


def listdir_fp(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


def get_cache(args):
    """First, checks to see if args.cache is specified. If not, it will then
    look for the HDSPIN_CACHE_DIR environment variable. If that also does not
    exist, then it will raise a RuntimeError. Returns the directory location.
    """

    if args.cache is not None:
        return args.cache

    env_var = os.environ.get("HDSPIN_CACHE_DIR", None)
    if env_var is not None:
        return env_var

    raise RuntimeError("Unknown cache location.")


def make_basename(nspin, beta, bc, dynamics, landscape, timesteps, loopN):
    """Creates the base filename for the simulation by simply combining the
    number of spins, beta, beta critical, the dynamics flag, the landscape
    flag, and the number of timesteps. This is a simple helper to ensure
    consistency."""

    return f"{nspin}_{beta:.03f}_{bc:.03f}_{dynamics}_{landscape}_{timesteps}_{loopN}"


def make_grids(args, grid_path):
    """Writes the grids to disk in the standard spot."""

    config = yaml.safe_load(open("configs/grid_configs/config.yaml"))

    nMC = int(10**args.timesteps)

    # Quick check to ensure that the number of nMC steps fits into a long long
    assert nMC < MAX_LONG

    energy_grid = np.unique(np.logspace(
        0, np.log10(nMC), config['energy_gridpoints'],
        dtype=int, endpoint=True
    ))

    tw_max = nMC // (args.dw + 1.0)

    # Define the first grid.
    pi_g1 = np.unique(np.logspace(
        0, np.log10(tw_max), config['pi_gridpoints'], dtype=int, endpoint=True
    ))

    # The second grid is directly related to the first via
    # tw -> tw + tw * dw
    pi_g2 = (pi_g1 * (args.dw + 1.0)).astype(int)

    np.savetxt(f"{grid_path}/energy.txt", energy_grid, fmt="%i")
    np.savetxt(f"{grid_path}/pi1.txt", pi_g1, fmt="%i")
    np.savetxt(f"{grid_path}/pi2.txt", pi_g2, fmt="%i")


def make_directory_and_configs(args):
    """Makes the appropriate directories for a single job. Also returns the
    directory name. Also returns the max index + 1 corresponding to the
    previously processed trials, so that we can actually continue to add
    statistics to previously run trials, if desired."""

    cache = args.cache
    basename = make_basename(
        args.nspin, args.beta, args.beta_critical, args.dynamics,
        args.landscape, args.timesteps, args.loopN
    )
    base_dir = os.path.join(cache, basename)

    max_index = 0
    if os.path.exists(base_dir):
        dirs_in_results = os.listdir(os.path.join(base_dir, "results"))
        dirs_in_results = [d for d in dirs_in_results if "energy" in d]
        if len(dirs_in_results) > 0:
            dirs_in_results = sorted(dirs_in_results)[-1]
            max_index = int(dirs_in_results.split("_energy.txt")[0]) + 1
            warnings.warn(
                f"Base directory {base_dir} exists; will resume "
                f"at {max_index} based on saved results like *_energy.txt",
                RuntimeWarning
            )

    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(os.path.join(base_dir, "scripts"), exist_ok=True)
    grid_path = os.path.join(base_dir, "grids")
    os.makedirs(grid_path, exist_ok=True)
    os.makedirs(os.path.join(base_dir, "results"), exist_ok=True)
    os.makedirs(os.path.join(base_dir, "final"), exist_ok=True)

    make_grids(args, grid_path)

    return base_dir, max_index


def write_bash_script(args, base_dir, max_index):
    """Writes a single bash script to the working directory, which stacks
    various primed protocols sequentially. This is intended for local debugging
    and should probably not be used for production runs."""

    # Ensure the total number of configs fits into a long long
    assert int(2**args.nspin) < MAX_LONG

    script_name = "local_submit.sh"
    results_dir = os.path.join(base_dir, "results")
    grids_dir = os.path.join(base_dir, "grids")
    timesteps = args.timesteps
    args_str = f"{results_dir} {grids_dir} {timesteps} {args.nspin} " \
        f"{args.beta} {args.beta_critical} {args.landscape} " \
        f"{args.dynamics} {int(args.loopN)} {max_index} 0 {args.nsim}"

    # Then we append
    if os.path.exists(script_name):
        with open(script_name, 'a') as f:
            f.write(f"./exe/main.out {args_str}\n")
        print(f"Script {script_name} appended with new trial")

    # Else we write a new file
    else:
        with open(script_name, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("\n")
            f.write(f"./exe/main.out {args_str}\n")
        print(f"Script {script_name} written to disk")


def write_SLURM_script(args, base_dir, max_index):
    """Writes the SLURM submission script to the cache directory pertaining
    to this run."""

    submit_fname = os.path.join(base_dir, "scripts/submit.sh")
    configs = yaml.safe_load(
        open(f"configs/slurm_configs/config.yaml")
    )

    partition = configs['partition']
    runtime = configs['time']
    account = configs['account']
    constraint = configs['constraint']
    n_cpus_per_job = 1
    total_cpus = configs['total_cpus']
    modules = configs['modules']
    max_concurrent = configs['max_concurrent']

    # simplify things: let's use an even divisor
    assert total_cpus % n_cpus_per_job == 0
    assert args.nsim % total_cpus == 0
    arr_len = total_cpus // n_cpus_per_job
    sims_per_job = args.nsim // arr_len

    results_dir = os.path.join(base_dir, "results")
    grids_dir = os.path.join(base_dir, "grids")
    timesteps = args.timesteps
    args_str = f"{results_dir} {grids_dir} {timesteps} {args.nspin} " \
        f"{args.beta} {args.beta_critical} {args.landscape} {args.dynamics} " \
        f"{int(args.loopN)} {max_index} $SLURM_ARRAY_TASK_ID {sims_per_job}"

    with open(submit_fname, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("\n")
        f.write(f"#SBATCH --job-name=hdspin \n")
        f.write(f"#SBATCH -p {partition}\n")
        if runtime is not None:
            f.write(f"#SBATCH -t {runtime}\n")
        if account is not None:
            f.write(f"#SBATCH --account={account}\n")
        if constraint is not None:
            f.write(f"#SBATCH -C {constraint}\n")

        f.write("#SBATCH -n 1\n")

        f.write(f"#SBATCH --cpus-per-task={n_cpus_per_job}\n")
        f.write(f"#SBATCH --mem-per-cpu={configs['mem_per_cpu']}M\n")

        f.write(f"#SBATCH --output=job_data/hdspin_%A_%a.out\n")
        f.write(f"#SBATCH --error=job_data/hdspin_%A_%a.err\n")

        if max_concurrent is None:
            f.write(f"#SBATCH --array=0-{arr_len - 1}\n")
        else:
            f.write(f"#SBATCH --array=0-{arr_len - 1}%{max_concurrent}\n")
        f.write('\n')

        if modules is not None:
            for m in modules:
                f.write(f"module load {m}\n")
            f.write('\n')

        f.write(f"export OMP_NUM_THREADS={n_cpus_per_job}\n")
        f.write('\n')

        f.write(f'./exe/main.out {args_str}\n')


def _call_subprocess(script):
    process = subprocess.Popen(
        script, shell=True, stdout=subprocess.PIPE, universal_newlines=True
    )
    process.wait()


def submit(args):
    """Submits all jobs in the cache to the job controller."""

    jobs = os.listdir(args.cache)
    print("Preparing to submit jobs:")
    for j in jobs:
        print(j)

    if not args.force:
        user_input = input("Continue? [yes]/no\n")
        if user_input != 'yes':
            print("Exiting")
            return

    print("Submitting jobs!")

    for j in jobs:
        script = os.path.join(args.cache, j, "scripts/submit.sh")
        _call_subprocess(f'mv {script} .')
        _call_subprocess("sbatch submit.sh")
        _call_subprocess(f'mv submit.sh {args.cache}/{j}/scripts')
