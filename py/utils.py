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

    resume_at = 0
    if os.path.exists(base_dir):
        dirs_in_results = os.listdir(os.path.join(base_dir, "results"))
        if len(dirs_in_results) > 0:
            dirs_in_results = sorted(dirs_in_results)[-1]
            resume_at = int(dirs_in_results.split(".txt")[0]) + 1
            raise RuntimeWarning(
                f"Base directory {base_dir} exists; will resume at {resume_at}"
            )

    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(os.path.join(base_dir, "scripts"), exist_ok=True)
    os.makedirs(os.path.join(base_dir, "results"), exist_ok=True)
    os.makedirs(os.path.join(base_dir, "final"), exist_ok=True)
    return base_dir, resume_at


def approximate_mem_per_cpu(args, cpu_per_task):
    """Using an empirical formula that the approximate required size of the
    arrays used during computation are ~8 (bytes) * 2 (arrays) * 2^N
    (configurations), we approximate the amount of memory allocated per CPU.
    Will return a minimum value of 2MB."""

    n_configs = int(2**args.nspin)
    approximate_memory = 32 * n_configs / 1e6  # MB

    # Add 10% overhead
    with_overhead = int(1.1 * approximate_memory / cpu_per_task)

    return max(10, with_overhead)


def write_SLURM_script(args, base_dir, resume_at):
    """Writes the SLURM submission script to the cache directory pertaining
    to this run."""

    submit_fname = os.path.join(base_dir, "scripts/submit.sh")
    configs = yaml.safe_load(open(f"slurm_configs/config.yaml"))

    partition = configs['partition']
    runtime = configs['time']
    account = configs['account']
    constraint = configs['constraint']
    n_cpu_per_job = configs['n_cpu_per_job']
    total_cpu = configs['total_cpu']
    module = configs['module']
    nsim = args.nsim

    # simplify things: let's use an even divisor
    assert total_cpu % n_cpu_per_job == 0
    assert nsim % total_cpu == 0
    arr_len = total_cpu // n_cpu_per_job
    sims_per_job = nsim // arr_len

    approximate_memory = approximate_mem_per_cpu(args, n_cpu_per_job)

    results_dir = os.path.join(base_dir, "results")
    timesteps = 10**args.timesteps
    args_str = f"{results_dir} {timesteps} {args.nspin} {args.beta} " \
        f"{args.beta_critical} {args.landscape} {args.dynamics} {resume_at} " \
        f"$SLURM_ARRAY_TASK_ID {sims_per_job}"

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
        f.write(f"#SBATCH -c {n_cpu_per_job}\n")
        f.write(f"#SBATCH --mem-per-cpu={approximate_memory}MB\n")

        f.write(f"#SBATCH --output=job_data/hdspin_%A_%a.out\n")
        f.write(f"#SBATCH --error=job_data/hdspin_%A_%a.err\n")
        f.write('\n')

        f.write(f"#SBATCH --array=0-{arr_len - 1}\n")
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
