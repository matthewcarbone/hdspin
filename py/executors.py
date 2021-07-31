__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"


import copy
from datetime import datetime
from itertools import product
from math import floor
import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil
import yaml

from py import utils as u


MAX_LONG = 9223372036854775807


def merge_dictionaries(dict1, dict2):
    return {**dict1, **dict2}


def execution_parameters_permutations(dictionary):
    """Inputs a dictionary of a format such as

    eg = {
        hp1: [1, 2]
        hp2: [3, 4]
    }

    and returns a list of all permutations:

    eg1 = {
        hp1: 1
        hp2: 3
    }

    eg2 = {
        hp1: 1
        hp2: 4
    }

    eg3 = {
        hp1: 2
        hp2: 3
    }

    eg4 = {
        hp1: 2
        hp2: 4
    }
    """

    return [
        dict(zip(dictionary, prod)) for prod in product(
            *(dictionary[ii] for ii in dictionary)
        )
    ]


def get_unique_dictionary_permuations(
    dictionary,
    key_to_permute='permute_over',
    key_no_permute='no_permute_over'
):
    """Returns a list of dictionaries, permuting every combination in
    key_to_permute."""

    perms = execution_parameters_permutations(dictionary[key_to_permute])
    other = dictionary[key_no_permute]
    dictionary_copied = copy.deepcopy(dictionary)
    tmp = [merge_dictionaries(other, perm) for perm in perms]
    del dictionary_copied[key_to_permute]
    del dictionary_copied[key_no_permute]
    return [merge_dictionaries(tt, dictionary_copied) for tt in tmp]


def cleanup():
    """Removes all run files, saved data in the cache and auxiliary files
    in the local directory."""

    p = Path("local_submit.sh")
    if p.exists():
        p.unlink()
        print(f"Removed: {str(p)}")

    p = Path("main.out")
    if p.exists():
        p.unlink()
        print(f"Removed: {str(p)}")

    p = Path(".FIFO.yaml")
    if p.exists():
        p.unlink()
        print(f"Removed: {str(p)}")

    p = Path("job_data")
    if p.is_dir():
        shutil.rmtree(str(p))
        print(f"Removed tree: {str(p)}")

    p = u.get_cache()
    if p is not None:
        p = Path(p)
    if p.exists():
        shutil.rmtree(p)
        print(f"Removed tree: {str(p)}")

    print("Done")


class SingleRankParameters:
    """Parameter container for a single run.

    Parameters
    ----------
    _tmp_data : dict
    """

    @staticmethod
    def _assert_positive_integer(ii):
        assert isinstance(ii, int)
        assert ii > 0

    def _set_log10_timesteps(self, _tmp_data):
        self.log10_timesteps = _tmp_data['log10_timesteps']
        SingleRankParameters._assert_positive_integer(self.log10_timesteps)

    def _set_n_tracers(self, _tmp_data):
        self.n_tracers = _tmp_data['n_tracers']
        SingleRankParameters._assert_positive_integer(self.n_tracers)

    def _set_n_spins(self, _tmp_data):
        self.n_spins = _tmp_data['n_spins']
        SingleRankParameters._assert_positive_integer(self.n_spins)

    def _set_dynamics(self, _tmp_data):
        self.dynamics = _tmp_data['dynamics']
        assert isinstance(self.dynamics, str)
        assert self.dynamics in [
            'standard', 'standard-loop', 'standard-divN', 'gillespie',
            'gillespie-divN'
        ]

    def _set_landscape(self, _tmp_data):
        self.landscape = _tmp_data['landscape']
        assert isinstance(self.landscape, str)
        assert self.landscape in ['erem', 'rem']

    def _set_beta(self, _tmp_data):
        beta = _tmp_data.get('beta')
        temperature = _tmp_data.get('temperature')
        if beta is None:
            assert temperature is not None
            assert temperature > 0.0
            self.beta = 1.0 / temperature
        elif temperature is None:
            assert beta is not None
            assert beta > 0.0
            self.beta = beta
        else:
            raise RuntimeError("Must set beta OR temperature, not both")

    def _set_critical_beta(self, _tmp_data):
        """Note this might depend on the dynamics, which should be set
        previously."""

        beta_crit = _tmp_data['beta_c']
        if beta_crit is not None:
            self.beta_c = beta_crit
            assert self.beta_c > 0.0
            return

        if self.landscape == 'erem':
            self.beta_c = 1.0
        elif self.landscape == 'rem':

            # This is ~sqrt(2 ln 2)
            self.beta_c = 1.177410022515475
        else:
            raise RuntimeError("Unknown landscape in setting beta_critical")

    def _set_aging_dw(self, _tmp_data):
        self.aging_dw = _tmp_data['aging_dw']
        if self.aging_dw is None:
            self.aging_dw = 0.5
            return

        assert self.aging_dw > 0.0

    def __init__(self, _tmp_data=None, path=None):

        if _tmp_data is None:
            assert path is not None
            _tmp_data = yaml.safe_load(open(path, 'r'))

        # Set the appropriate class attributes after running the necessary
        # sanity checks. We start with the runtime parameters:
        self._set_log10_timesteps(_tmp_data)
        self._set_n_tracers(_tmp_data)
        self._set_n_spins(_tmp_data)
        self._set_dynamics(_tmp_data)
        self._set_landscape(_tmp_data)
        self._set_beta(_tmp_data)
        self._set_critical_beta(_tmp_data)
        self._set_aging_dw(_tmp_data)

    def filename(self):
        """Returns the filename for the simulation by simply combining the
        number of spins, beta, beta critical, the dynamics flag, the landscape
        flag, and the number of timesteps. This is a simple helper to ensure
        consistency."""

        return f"{self.dynamics}_{self.landscape}_{self.log10_timesteps}" + \
            f"_{self.n_spins}_{int(self.beta * 1000)}"

    def arg_string(self, results_dir, grids_dir, final_dir):
        """Returns the arguments needed to feed to the C++ MPI code."""

        return f"{results_dir} {grids_dir} {final_dir} " \
            f"{self.log10_timesteps} {self.n_spins} {self.beta} " \
            f"{self.beta_c} {self.landscape} {self.dynamics} {self.n_tracers}"

    def save(self, path):
        """Saves a yaml file of these parameters to disk at the specified
        path."""

        with open(path, 'w') as f:
            yaml.safe_dump(self.__dict__, f)


class Primer:
    """
    Parameters
    ----------
    config_location : str
        The absolute path to the yaml file containing the input parameters.
        It corresponds to possibly multiple MPI runs.

    Attributes
    ----------
    all_params : list
        A list of SingleRankParameters classes containing all of the
        runtime (and other) parameters.
    """

    def _set_cache(self):
        self.cache = u.get_cache()
        if not Path(self.cache).exists():
            Path(self.cache).mkdir(parents=True, exist_ok=True)
            assert Path(self.cache).exists()
        Path("job_data").mkdir(exist_ok=True)

    def _set_grid_info(self, _tmp_data):
        self.energy_gridpoints = _tmp_data['energy_gridpoints']
        SingleRankParameters._assert_positive_integer(self.energy_gridpoints)
        self.pi_gridpoints = _tmp_data['pi_gridpoints']
        SingleRankParameters._assert_positive_integer(self.pi_gridpoints)

    def __init__(self, config_location, fifo_queue_path=".FIFO.yaml"):
        self.fifo_queue_path = Path(fifo_queue_path)
        _params = yaml.safe_load(open(config_location, 'r'))
        _uniqe_params = get_unique_dictionary_permuations(_params)
        self.all_params = [SingleRankParameters(xx) for xx in _uniqe_params]
        self._set_cache()
        self._set_grid_info(_params)
        print(f"\nPreparing {len(self.all_params)} MPI trial(s)")
        print(f"Cache directory: {self.cache}")
        print("Summary (note not all parameters are displayed):\n")
        p = pd.DataFrame([c.__dict__ for c in self.all_params])
        print(p, '\n')

    def make_grids(self, dw, log10_timesteps, grid_path):
        """Writes the grids to disk in the standard spot."""

        nMC = int(10**log10_timesteps)

        # Quick check to ensure that the number of nMC steps fits into a long
        # long
        assert nMC < MAX_LONG

        energy_grid = np.unique(np.logspace(
            0, np.log10(nMC), self.energy_gridpoints, dtype=int, endpoint=True
        ))

        # Adds 0 to the energy grid so we can record the value of the energy
        # at 0 itself.
        energy_grid = np.array([0] + list(energy_grid))

        tw_max = nMC // (dw + 1.0)

        # Define the first grid.
        pi_g1 = np.unique(np.logspace(
            0, np.log10(tw_max), self.pi_gridpoints, dtype=int, endpoint=True
        ))

        # The second grid is directly related to the first via
        # tw -> tw + tw * dw
        pi_g2 = (pi_g1 * (dw + 1.0)).astype(int)

        np.savetxt(f"{grid_path}/energy.txt", energy_grid, fmt="%i")
        np.savetxt(f"{grid_path}/pi1.txt", pi_g1, fmt="%i")
        np.savetxt(f"{grid_path}/pi2.txt", pi_g2, fmt="%i")

    def _append_queue(self, path):
        """For convenience, so the user doesn't need to copy/paste the package
        names each time they submit a job. This is a FIFO (last in first out)
        queue, so the most recent primed job will be the one submitted if the
        user does not specify a package."""

        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

        if self.fifo_queue_path.exists():
            queue = yaml.safe_load(open(str(self.fifo_queue_path)))
            queue['trials'].append({
                'date_primed': dt_string,
                'path': str(path)
            })
        else:
            queue = {
                'trials': [{
                    'date_primed': dt_string,
                    'path': str(path)
                }]
            }

        with open(str(self.fifo_queue_path), 'w') as f:
            yaml.dump(queue, f, default_flow_style=False)

    def prime(self):
        """Makes a directory for each trial."""

        now = datetime.now()
        dt_string = now.strftime("%Y%m%d%H%M%S")

        for srp in self.all_params:  # Single Rank Parameter (srp)
            fname = Path(srp.filename() + f"_{dt_string}")
            fname = Path(self.cache) / fname
            self._append_queue(fname)
            fname.mkdir(parents=True, exist_ok=False)
            grid_path = fname / Path("grids")
            grid_path.mkdir(parents=False, exist_ok=False)
            self.make_grids(srp.aging_dw, srp.log10_timesteps, grid_path)
            (fname / Path("final")).mkdir(parents=False, exist_ok=False)
            (fname / Path("results")).mkdir(parents=False, exist_ok=False)
            srp.save(str(fname / Path("params.yaml")))
            print(f"Made directories: {str(fname)}")


class Executor:
    """Runs the code!"""

    def _set_cache(self):
        self.cache = u.get_cache()
        assert Path(self.cache).exists()

    def _load_from_FIFO_queue(self, load_all=False):
        """Loads the last entry in the saved FIFO queue, returns that entry,
        and pops that entry from the list, resaving the file."""

        queue = yaml.safe_load(open(self.fifo_queue_path))

        if len(queue['trials']) < 1:
            print("FIFO queue is empty, run prime before execute. Exiting.")
            self.fifo_queue_path.unlink()
            exit(0)

        if not load_all:

            # Get the first entry
            last_entry = queue['trials'][0]

            # Resave the file
            queue['trials'] = queue['trials'][1:]
            with open(self.queue_path, 'w') as f:
                yaml.dump(queue, f, default_flow_style=False)

            return [last_entry['path']]

        self.fifo_queue_path.unlink()
        return [entry['path'] for entry in queue['trials']]

    def write_SLURM_script(self, base_dir, srp):
        """Writes the SLURM submission script to the cache directory pertaining
        to this run."""

        base_dir = Path(base_dir)

        partition = self.slurm_params.get("partition")
        runtime = self.slurm_params.get("time")
        account = self.slurm_params.get("account")
        constraint = self.slurm_params.get("constraint")
        N_nodes = self.slurm_params['N_nodes']
        mem = self.slurm_params.get("mem_per_node")
        SBATCH_lines = [
            "#!/bin/bash\n\n",
            "#SBATCH --job-name=hdspin\n",
            f"#SBATCH -p {partition}\n" if partition is not None else "",
            f"#SBATCH -t {runtime}\n" if runtime is not None else "",
            f"#SBATCH --account={account}\n" if account is not None else "",
            f"#SBATCH -C {constraint}\n" if constraint is not None else "",
            f"#SBATCH -N {N_nodes}\n",
            f"#SBATCH --mem={mem}G\n" if mem is not None else "",
            "\n",
            "#SBATCH --output=job_data/hdspin_%A_%a.out\n",
            "#SBATCH --error=job_data/hdspin_%A_%a.err\n",
            "\n"
        ]

        # Parse the MPI information
        cores_per_node = self.slurm_params["cores_per_node"]
        hyperthreads_per_core = self.slurm_params["hyperthreads_per_core"]
        tasks_per_node = self.slurm_params["tasks_per_node"]
        assert isinstance(tasks_per_node, int)
        SBATCH_lines.append(f"#SBATCH --tasks-per-node={tasks_per_node}\n")
        phys_cores_per_task = int(floor(cores_per_node/tasks_per_node))
        c = phys_cores_per_task * hyperthreads_per_core
        SBATCH_lines.append(f"#SBATCH -c {c}\n\n")

        # Handle the modules
        for line in self.slurm_params['other_lines']:
            SBATCH_lines.append(f"{line}\n")
        SBATCH_lines.append("\n")

        # Handle threads optionally
        max_threads = self.slurm_params.get("max_threads")
        if max_threads is None:
            threads = phys_cores_per_task
        else:
            threads = min(phys_cores_per_task, max_threads)
        SBATCH_lines.append(f"export OMP_NUM_THREADS={threads}\n")

        results_dir = str(base_dir / Path("results"))
        grids_dir = str(base_dir / Path("grids"))
        final_dir = str(base_dir / Path("final"))
        args_str = srp.arg_string(results_dir, grids_dir, final_dir)

        submit_fname = base_dir / Path("submit.sh")
        with open(submit_fname, 'w') as f:
            for line in SBATCH_lines:
                f.write(line)
            f.write(f'\nmpiexec ./exe/main.out {args_str}\n')
            # f.write(
            #     f"python3 run.py eval --directory {base_dir}\n"
            # )

    def write_bash_script(args, base_dir, srp):
        """Writes a single bash script to the working directory, which stacks
        various primed protocols sequentially. This is intended for local
        debugging and should probably not be used for production runs."""

        results_dir = str(Path(base_dir) / Path("results"))
        grids_dir = str(Path(base_dir) / Path("grids"))
        final_dir = str(base_dir / Path("final"))
        args_str = srp.arg_string(results_dir, grids_dir, final_dir)

        script_name = "local_submit.sh"

        # Then we append
        if os.path.exists(script_name):
            with open(script_name, 'a') as f:
                f.write(f"mpiexec -np 6 ./exe/main.out {args_str}\n")
            print(f"Script {script_name} appended with new trial")

        # Else we write a new file
        else:
            with open(script_name, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("\n")
                f.write(f"mpiexec -np 6 ./exe/main.out {args_str}\n")
            print(f"Script {script_name} written to disk")

    def __init__(
        self, config_location, run_one, local=False,
        fifo_queue_path=".FIFO.yaml", exe_path=Path("exe/main.out")
    ):
        if not exe_path.exists():
            raise RuntimeError("Run Make before using run.py execute")
        self._set_cache()
        self.fifo_queue_path = Path(fifo_queue_path)
        self.slurm_params = yaml.safe_load(open(config_location, 'r'))['SLURM']
        self.config_location = Path(config_location)
        self.primed = self._load_from_FIFO_queue(load_all=not run_one)
        for base_dir in self.primed:
            p = Path(base_dir) / Path("params.yaml")
            params = yaml.safe_load(open(str(p), 'r'))
            srp = SingleRankParameters(params)
            if local:
                self.write_bash_script(base_dir, srp)
            else:
                self.write_SLURM_script(base_dir, srp)
        print(f"{len(self.primed)} jobs with local={local}: READY")

    def execute(self):
        """Submits the SLURM jobs."""

        cache_full_paths = u.listdir_fp(self.cache)
        for path in cache_full_paths:
            script_loc = Path(path) / Path("submit.sh")
            u.run_command(f"mv {script_loc} .", silent=False)
            u.run_command("sbatch submit.sh", silent=False)
            u.run_command(f"mv submit.sh {path}", silent=False)
