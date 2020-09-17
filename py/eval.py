#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import os

import numpy as np
import pickle

from py import utils as u


class Evaluator:
    """Evaluates all available results.

    Parameters
    ----------
    args
        As passed via the argparser.

    Attributes
    ----------
    cache : str
        The location to the cache directory containing all trials.
    all_dirs : list
        A list of the full paths to directories like 16_0.750_1.000_0_0_6,
        containing results for a run corresponding to the parameters as
        specified in the filename.
    """

    def __init__(self, args):
        """Initializer."""

        self.cache = u.get_cache(args)
        _all_dirs = u.listdir_fp(self.cache)
        self.all_dirs = [d for d in _all_dirs if os.path.isdir(d)]

    def eval_aging_config(self):
        """Evalutes all saved pi/aging config results."""

        print(
            f"Pi/Aging Config: evaluating {len(self.all_dirs)} total "
            "directories"
        )

        for full_dir_path in self.all_dirs:

            results_path = os.path.join(full_dir_path, 'results')
            all_trials = os.listdir(results_path)
            all_trials.sort()
            pi1 = np.array([
                np.loadtxt(os.path.join(results_path, f), delimiter=" ")
                for f in all_trials if "pi1_config" in f
            ])
            pi2 = np.array([
                np.loadtxt(os.path.join(results_path, f), delimiter=" ")
                for f in all_trials if "pi2_config" in f
            ])

            # Find where the config index, the standard config integer rep, and
            # the inherent structure rep are equal
            eq = pi1 == pi2

            final_path = os.path.join(
                full_dir_path, 'final/aging_config_sd.pkl'
            )

            print(f"Evaluating {full_dir_path} -> {final_path}")

            pickle.dump(eq, open(final_path, 'wb'), protocol=4)

        print("Done\n")

    def eval_psi_config(self):
        """Evaluates all saved psi config results."""

        print(
            f"Psi Config: evaluating {len(self.all_dirs)} total directories"
        )

        for full_dir_path in self.all_dirs:

            results_path = os.path.join(full_dir_path, 'results')
            all_trials = os.listdir(results_path)
            res = [
                np.loadtxt(os.path.join(results_path, f), delimiter=" ")
                for f in all_trials if "psi_config" in f
            ]

            # Construct dictionaries out of each result while also keeping
            # track of the max key value. There's probably a way to do this
            # using list comprehension, but this is easier for now, and the
            # total number of results will never more more than ~10k, so it
            # should be ok.
            dict_res = []
            dict_res_IS = []
            max_key = 0
            for res_arr in res:

                # Standard
                d = {
                    int(key): int(value) for key, value
                    in zip(res_arr[:, 0], res_arr[:, 1])
                }
                for key, value in d.items():
                    if key > max_key and value > 0:
                        max_key = key
                dict_res.append(d)

                # Inherent structure
                d = {
                    int(key): int(value) for key, value
                    in zip(res_arr[:, 0], res_arr[:, 2])
                }
                for key, value in d.items():
                    if key > max_key and value > 0:
                        max_key = key
                dict_res_IS.append(d)

            # Standard
            stats = np.zeros(shape=(len(dict_res), max_key + 1))
            for ii, d in enumerate(dict_res):
                for key, value in d.items():
                    if value > 0:
                        stats[ii, key] = value
            final_path = os.path.join(full_dir_path, 'final/psi_config.txt')
            print(f"Evaluating {full_dir_path} - shape {stats.shape}")
            np.savetxt(final_path, stats)

            # Inherent structure
            stats_IS = np.zeros(shape=(len(dict_res_IS), max_key + 1))
            for ii, d in enumerate(dict_res_IS):
                for key, value in d.items():
                    if value > 0:
                        stats_IS[ii, key] = value
            final_path = os.path.join(full_dir_path, 'final/psi_config_IS.txt')
            print(f"Evaluating {full_dir_path} -> {final_path}")
            np.savetxt(final_path, stats_IS)

        print("Done\n")

    def eval_traj(self):
        """Performs the evaluate of the energy, trajectories and inherent
        structure for every directory as listed in all_dirs."""

        print(f"Energy: evaluating {len(self.all_dirs)} total directories")

        for full_dir_path in self.all_dirs:

            results_path = os.path.join(full_dir_path, 'results')
            all_trials = os.listdir(results_path)
            res = np.array([
                np.loadtxt(os.path.join(results_path, f), delimiter=" ")
                for f in all_trials if "energy" in f
            ])
            # print(f"Read trials of shape {res.shape} from {results_path}")

            # Assert that the time-grids are identical for all trials
            assert len(np.unique(res[:, :, 0])) == res.shape[1]

            # Load in the energy and inherent structure energy
            energies = res[:, :, 2]
            energies_inherent = res[:, :, 4]

            # print("\n\n")
            # print(energies.min())
            # print(energies_inherent.min())
            # print("\n\n")

            # print(f"Loaded energies of shape {energies.shape}")

            final_path = os.path.join(full_dir_path, 'final/energy.txt')

            to_save = np.array([
                res[0, :, 0].squeeze(),
                energies.mean(axis=0).squeeze(),
                energies.std(axis=0).squeeze(),
                energies_inherent.mean(axis=0).squeeze(),
                energies_inherent.std(axis=0).squeeze()
            ])

            print(f"Evaluating {full_dir_path} -> {final_path}")

            np.savetxt(final_path, to_save.T)

        print("Done\n")
