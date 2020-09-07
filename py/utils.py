#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"


import os
import subprocess
import yaml


def make_basename(nspin, beta, bc, dynamics, landscape, timesteps):
    """Creates the base filename for the simulation by simply combining the
    number of spins, beta, beta critical, the dynamics flag, the landscape
    flag, and the number of timesteps. This is a simple helper to ensure
    consistency."""

    return f"{nspin}_{beta:.03f}_{bc:.03f}_{dynamics}_{landscape}_{timesteps}"


def make_directory(args):
    """Makes the appropriate directories for a single job. Also returns the
    directory name. Also returns the max index + 1 corresponding to the
    previously processed trials, so that we can actually continue to add
    statistics to previously run trials, if desired."""

    cache = args.cache
    basename = make_basename(
        args.nspin, args.beta, args.beta_critical, args.dynamics,
        args.landscape, args.timesteps
    )
    base_dir = os.path.join(cache, basename)

    max_index = 0
    if os.path.exists(base_dir):
        dirs_in_results = os.listdir(os.path.join(base_dir, "results"))
        if len(dirs_in_results) > 0:
            dirs_in_results = sorted(dirs_in_results)[-1]
            max_index = int(dirs_in_results.split(".txt")[0]) + 1
            raise RuntimeWarning(
                f"Base directory {base_dir} exists; will resume at {max_index}"
            )

    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(os.path.join(base_dir, "scripts"), exist_ok=True)
    os.makedirs(os.path.join(base_dir, "results"), exist_ok=True)
    os.makedirs(os.path.join(base_dir, "final"), exist_ok=True)
    return base_dir, max_index


def approximate_mem_per_cpu(args, cpu_per_task):
    """Using an empirical formula that the approximate required size of the
    arrays used during computation are ~8 (bytes) * 2 (arrays) * 2^N
    (configurations), we approximate the amount of memory allocated per CPU.
    Will return a minimum value of 2MB."""

    n_configs = int(2**args.nspin)
    approximate_memory = 16 * n_configs / 1e6  # MB

    # Add 5% overhead
    with_overhead = int(1.05 * approximate_memory)

    # Divide by the number of cpus per MPI task
    mem = with_overhead // cpu_per_task
    return max(2, mem)


def write_SLURM_script(args, base_dir, max_index):
    """Writes the SLURM submission script to the cache directory pertaining
    to this run."""

    submit_fname = os.path.join(base_dir, "scripts/submit.sh")
    configs = yaml.safe_load(
        open(f"slurm_configs/{args.slurm_cluster_preset}.yaml")
    )

    partition = configs['partition']
    runtime = configs['time']
    account = configs['account']
    n_tasks = configs['n_tasks']
    n_cpus_per_task = configs['n_cpus_per_task']
    module = configs['module']

    approximate_memory = approximate_mem_per_cpu(args, n_cpus_per_task)

    results_dir = os.path.join(base_dir, "results")
    timesteps = 10**args.timesteps
    args_str = f"{results_dir} {timesteps} {args.nspin} {args.beta} " \
        f"{args.beta_critical} {args.landscape} {args.dynamics} {args.nsim} " \
        f"{max_index}"

    with open(submit_fname, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("\n")
        f.write(f"#SBATCH --job-name=hdspin \n")
        f.write(f"#SBATCH -p {partition}\n")
        if runtime is not None:
            f.write(f"#SBATCH -t {runtime}\n")
        if account is not None:
            f.write(f"#SBATCH --account={account}\n")

        f.write(f"#SBATCH --ntasks={n_tasks}\n")
        f.write(f"#SBATCH --cpus-per-task={n_cpus_per_task}\n")
        f.write(f"#SBATCH --mem-per-cpu={approximate_memory}M\n")

        f.write(f"#SBATCH --output=job_data/hdspin_%A.out\n")
        f.write(f"#SBATCH --error=job_data/hdspin_%A.err\n")
        f.write('\n')

        if module is not None:
            f.write(f"module load {module}\n")
            f.write('\n')

        f.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
        f.write('\n')

        f.write(f'mpirun ./main.out {args_str}\n')


def _call_subprocess(script):
    process = subprocess.Popen(
        script, shell=True, stdout=subprocess.PIPE, universal_newlines=True
    )
    process.wait()


def submit(args):
    """Submits all jobs in the cache to the job controller."""

    # Move the script to the working directory
    script = os.path.join(args.cache, "scripts/submit.sh")
    _call_subprocess(f'mv {script} .')
    _call_subprocess("sbatch submit.sh")
    _call_subprocess(f'mv submit.sh {args.cache}/scripts')
