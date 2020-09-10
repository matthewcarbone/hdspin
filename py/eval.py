#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import os

import numpy as np

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
        # print(self.all_dirs)

    def eval_traj(self):
        """Performs the evaluate of the energy, trajectories and inherent
        structure for every directory as listed in all_dirs."""

        print(f"Evaluating {len(self.all_dirs)} total directories\n")

        for full_dir_path in self.all_dirs:
            print(f"Evaluating trajectories/energies for {full_dir_path}")

            results_path = os.path.join(full_dir_path, 'results')
            all_trials = os.listdir(results_path)
            res = np.array([
                np.loadtxt(os.path.join(results_path, f), delimiter=" ")
                for f in all_trials if "energy" in f
            ])
            print(f"Read trials of shape {res.shape} from {results_path}")

            # Assert that the time-grids are identical for all trials
            assert len(np.unique(res[:, :, 0])) == res.shape[1]

            # Load in the energy and inherent structure energy
            energies = res[:, :, 2].squeeze()
            energies_inherent = res[:, :, 4].squeeze()

            print(f"Loaded energies of shape {energies.shape}")

            final_path = os.path.join(full_dir_path, 'final/energy.txt')
            print(
                f"Saving energy and energy inherent structure to {final_path}"
            )

            to_save = np.array([
                res[0, :, 0].squeeze(), energies.mean(axis=0),
                energies.std(axis=0), energies_inherent.mean(axis=0),
                energies_inherent.std(axis=0)
            ])
            np.savetxt(final_path, to_save.T)
            print("Done\n")
