#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import yaml
import argparse
import subprocess
import os
import pickle
import numpy as np
from itertools import product

PARAMS_DIR = 'params/'
FINAL_RESULTS_DIR = 'final_results/'
DIRS_TO_MAKE = [
    'energy',
    'max_energy',
    'psi_c',
    'psi_b_S',
    'psi_b_E',
    'pi_b_S_1',
    'pi_b_S_2',
    'pi_b_E_1',
    'pi_b_E_2',
    'pi_in_b_S_1',
    'pi_in_b_S_2',
    'pi_in_b_E_1',
    'pi_in_b_E_2',
]


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

    combinations = [dict(zip(dictionary, prod))
                    for prod in product(*(dictionary[ii]
                                          for ii in dictionary))]
    return combinations


def merge_dict(dict1, dict2):
    return {**dict1, **dict2}


def parser():
    """Uses argparse to parse command line arguments and returns said
    arguments."""

    ap = argparse.ArgumentParser()

    ap.add_argument('--dryrun', action='store_true', dest='dryrun',
                    default=False, help='run in dryrun mode')

    ap.add_argument('--local', action='store_true', dest='local',
                    default=False, help='run locally')

    ap.add_argument('--concat', action='store_true', dest='concat',
                    default=False, help='concatenate the results')

    return ap.parse_args()


def make_directory_paths(list_of_parameters):
    """Creates the necessary directories if they do not exist."""

    for d1 in list_of_parameters:
        for d2 in DIRS_TO_MAKE:
            os.makedirs(
                os.path.join(d1['file_dump_loc'], d2),
                exist_ok=True
            )


def get_all_abs_paths(in_dir):
    """Gets all the absolute paths of the files in the argument."""

    _l = []
    for _dir in os.listdir(in_dir):
        _l.append(os.path.join(PARAMS_DIR, _dir))
    return _l


def print_dictionary_no_depth(d):
    """Prints the contents of a dictionary, `d` at a single level of depth."""

    for key, value in d.items():
        print(f"\t ~ {key} -> {value}")


def process_parameter_list(list_of_parameters):
    """The parameter input files cycle over a list of n_spins and beta values.
    This function processes parameter files of this form and returns the true
    list of the parameter files to run."""

    final_list = []

    for p in list_of_parameters:
        no_cycle_params = p['no_cycle']
        cycle_params = execution_parameters_permutations(p['cycle'])

        for cyc in cycle_params:
            new_dict = merge_dict(no_cycle_params, cyc)
            f = new_dict['file_dump_loc']
            dy = new_dict['dynamics']
            N = new_dict['n_spins']
            beta = new_dict['beta']
            new_dict['file_dump_loc'] = f"{f}/{dy}_{N:02}_{beta:.02f}"
            new_dict['final_name'] = f"{dy}_{N:02}_{beta:.02f}.pkl"
            final_list.append(new_dict)

    return final_list


def submit(list_of_parameters, local=False, dryrun=True):
    """Takes a list of parameter files and submits the jobs."""

    for cc, p in enumerate(list_of_parameters):

        # Get the dynamics switch
        if p['dynamics'] == 'standard':
            p['dynamics'] = 0
        elif p['dynamics'] == 'gillespie':
            p['dynamics'] = 1
        else:
            raise RuntimeError(f"Unknown dynamics switch {p['dynamics']}")

        # Get the landscape
        if p['landscape'] == 'erem':
            p['landscape'] = 0
            p['beta_critical'] = 1.0
        elif p['landscape'] == 'rem':
            p['landscape'] = 1
            p['beta_critical'] = np.sqrt(2.0 * np.log(2.0))
        else:
            raise RuntimeError(f"Unknown landscape switch {p['landscape']}")

        # Handle the scaling factor if desired
        if p['log10nMC'] == -1:
            exp_arg = p['beta_critical'] * p['beta'] * p['n_spins']
            p['log10nMC'] = np.log10(10.0 * np.exp(exp_arg))

        print(f"Submitting {cc:03}")
        print_dictionary_no_depth(p)
        print("\n")

        if dryrun:
            continue

        for ii in range(p['n_cpu']):

            subprocess.run([
                "sbatch" if not local else "bash",
                "_submit.sh",
                str(p['dynamics']),
                str(p['file_dump_loc']),
                str(p['n_spins']),
                str(p['beta']),
                str(p['beta_critical']),
                str(p['log10nMC']),
                str(p['landscape']),
                str(p['n_serial_runs']),
                str(ii)  # CPU index
            ])

            if local:
                break

        if local:
            break


def concat_results(directory, check_binary=False):
    """Concats the energy results."""

    files_in_dir = os.listdir(directory)
    res = []
    for f in files_in_dir:
        x = np.loadtxt(open(os.path.join(directory, f)))
        if check_binary:
            if np.any((x != 1) & (x != 0)):
                print(f"Warning: file {f} != 0|1")
        res.append(x)
    return np.array(res)


def process_pi_trackers(
    index_directory_1,
    index_directory_2
):
    """The index_directory points to the location of the files that index which
    basin the tracer is in at the recorded time, and the in_basin_directory
    indexes whether or not the tracer is in a basin at that time."""

    index1 = concat_results(index_directory_1)
    index2 = concat_results(index_directory_2)

    # Get the locations of the indexes where the tracer is in a basin at
    # the first time point:
    # same_and_in_basin = (index1 == index2) & (in1 == in2) & (in2 == 1)

    # Get the normalization, which is given by the number of tracers in a basin
    # at the first time
    # norm = (in1 == 1).sum(axis=0)

    # At both times the tracer must be in the *same* basin
    # cond1 = index1 == index2

    # At both times the tracer must be in a basin
    # cond2 = (in1 == in2) & (in2 == 1)

    # The union of both of these conditions is the result
    # overall_condition = cond1 & cond2
    return index1, index2


def concat(list_of_parameters):
    """Takes the results for each of the runs and amalgamates them along with
    errors. Saves results as a single pickle file in the relevant
    directory. The loaded results will always have the shape of
    (n_trials, n_sample) and should be processed accordingly."""

    os.makedirs("final_results", exist_ok=True)

    for ii, d in enumerate(list_of_parameters):
        print(f"Concatenating results in {d['file_dump_loc']}")
        res_dict = {
            'grids': dict(),
            'results': dict()
        }

        # Energy --------------------------------------------------------------
        res_dict['grids']['energy'] = np.loadtxt(
            open(os.path.join(d['file_dump_loc'], 'energy_grid.txt'))
        )
        res_dict['results']['energy'] = concat_results(
            os.path.join(d['file_dump_loc'], 'energy')
        )

        # Psi -----------------------------------------------------------------
        res_dict['results']['psi_config'] = concat_results(
            os.path.join(d['file_dump_loc'], 'psi_c')
        )

        # The Psi grid is a simple base2 logarithmic scale
        len_psi_grid = res_dict['results']['psi_config'].shape[1]
        res_dict['grids']['psi'] = [2**n for n in range(len_psi_grid)]

        res_dict['results']['psi_basin_S'] = concat_results(
            os.path.join(d['file_dump_loc'], 'psi_b_S')
        )
        res_dict['results']['psi_basin_E'] = concat_results(
            os.path.join(d['file_dump_loc'], 'psi_b_E')
        )

        # Pi ------------------------------------------------------------------
        res_dict['grids']['pi_1'] = np.loadtxt(
            open(os.path.join(d['file_dump_loc'], 'Pi_basin_grid_1.txt'))
        )
        res_dict['grids']['pi_2'] = np.loadtxt(
            open(os.path.join(d['file_dump_loc'], 'Pi_basin_grid_2.txt'))
        )
        res = process_pi_trackers(
            os.path.join(d['file_dump_loc'], 'pi_b_S_1'),
            os.path.join(d['file_dump_loc'], 'pi_b_S_2')
        )
        res_dict['results']['Pi_basin_S'] = res
        res = process_pi_trackers(
            os.path.join(d['file_dump_loc'], 'pi_b_E_1'),
            os.path.join(d['file_dump_loc'], 'pi_b_E_2')
        )

        #res_dict['grids']['energy_div_2'] = np.loadtxt(
        #    open(os.path.join(d['file_dump_loc'], 'energy_grid_t_div_2.txt'))
        #)
        #res_dict['results']['max_energy'] = concat_results(
        #    os.path.join(d['file_dump_loc'], 'max_energy')
        #)

        res_dict['results']['Pi_basin_E'] = res
        res_dict['parameters'] = d

        # Save results
        path = os.path.join(FINAL_RESULTS_DIR, d['final_name'])
        pickle.dump(res_dict, open(path, 'wb'), protocol=4)

        continue


# Main method -----------------------------------------------------------------

if __name__ == '__main__':

    # Parse the arguments
    args = parser()

    # Get all the parameter files
    all_paths = get_all_abs_paths(PARAMS_DIR)

    # Load in all the parameter files
    list_of_parameters = [yaml.safe_load(open(f, 'r')) for f in all_paths]
    list_of_parameters = process_parameter_list(list_of_parameters)

    if args.concat:
        concat(list_of_parameters)
        exit(0)

    # Check that the executable exists
    if not os.path.exists("main.o"):
        raise RuntimeError("Run `make` before submitting!")

    print("\nParameters to submit are:\n")
    for ii, d in enumerate(list_of_parameters):
        print(f"({ii + 1:03})")
        print_dictionary_no_depth(d)

    make_directory_paths(list_of_parameters)

    # Move scripts and run as appropriate
    subprocess.run(["mv", "scripts/_submit.sh", "."])

    if args.local:
        print("\nRunning only first parameter locally in serial...")
        submit(list_of_parameters, local=True, dryrun=args.dryrun)
    else:
        print("\nSubmitting to cluster...")
        submit(list_of_parameters, local=False, dryrun=args.dryrun)

    # Move the scripts back when done
    subprocess.run(["mv", "_submit.sh", "scripts/"])
