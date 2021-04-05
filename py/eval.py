#!/usr/bin/env python3

"""Contains all methods for evaluating the results of the hdspin calculations.
"""

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import os

import glob2
import numpy as np
from pathlib import Path

from py import utils as u


def process_file(file_name):
    """Reads the file_name using numpy.loadtxt with delimiter ' '.

    Parameters
    ----------
    file_name : str

    Returns
    -------
    np.ndarray
    """

    return np.loadtxt(file_name, delimiter=' ')


def read_files_via_numpy(file_list):
    """Reads all files as listed in the file_list using numpy.loadtxt. Can
    also use a multiprocessing Pool().map operation here but it runs into an
    OS error, indicating too many simultaneous open files.

    Parameters
    ----------
    file_list : list
        A list of full paths.

    Returns
    -------
    list
        A list of all the loaded numpy arrays.
    """

    return list(map(process_file, file_list))


class Evaluator:
    """Evaluates all available results.

    Parameters
    ----------
    args
        As passed via the argparser.

    Attributes
    ----------
    all_dirs : list
        A list of the full paths to directories like 16_0.750_1.000_0_0_6,
        containing results for a run corresponding to the parameters as
        specified in the filename.
    cache : str
        The location to the cache directory containing all trials.
    """

    @staticmethod
    def energy(root):
        """Performs the evaluate of the energy, trajectories and inherent
        structure for every directory as listed in all_dirs.

        Notes
        -----
        The loaded arrays have the following format
        COLUMN
            0: grid
            1: integer representation of the config
            2: energy
            3: integer representation of the inherent structure config
            4: inherent structure energy

        Parameters
        ----------
        root : str
            The base directory containing all of the results for this
            particular trial.
        """

        results_directory = Path(root) / Path("results")
        files = glob2.glob(str(results_directory / "*_energy.txt"))
        res = np.array(read_files_via_numpy(files))

        # Assert that the time-grids are identical for all trials
        assert len(np.unique(res[:, :, 0])) == res.shape[1]

        # Load in the energy and inherent structure energy. We don't really
        # care about the config indexes.
        energies = res[:, :, 2]
        energies_inherent = res[:, :, 4]
        N = res.shape[0]

        # Get the location of the directory to save the finalized results for
        # both the standard and...
        final_path = Path(root) / Path("final/energy.txt")
        sd = energies.std(axis=0).squeeze()
        to_save = np.array([
            res[0, :, 0].squeeze(),
            energies.mean(axis=0).squeeze(), sd, sd / np.sqrt(N - 1)
        ])
        np.savetxt(final_path, to_save.T)

        # ... inherent structure trajectories
        final_path = Path(root) / Path("final/energy_IS.txt")
        sd = energies_inherent.std(axis=0).squeeze()
        to_save = np.array([
            res[0, :, 0].squeeze(),
            energies_inherent.mean(axis=0).squeeze(), sd, sd / np.sqrt(N - 1)
        ])
        np.savetxt(final_path, to_save.T)

        # Also save some example raw trajectories.
        ntraj = min(res.shape[0], 50)
        final_path = Path(root) / Path("final/energy_eg_traj.txt")
        to_save = np.concatenate([
            res[0, :, 0][:, None], energies[:ntraj, :].T
        ], axis=1)
        np.savetxt(final_path, to_save)
        final_path = Path(root) / Path("final/energy_IS_eg_traj.txt")
        to_save = np.concatenate([
            res[0, :, 0][:, None], energies_inherent[:ntraj, :].T
        ], axis=1)
        np.savetxt(final_path, to_save)

        print(f"\tEnergy done")

    @staticmethod
    def psi_config(root, inherent_structure):
        """Evaluates all saved psi config results.

        Parameters
        ----------
        root : str
            The base directory containing all of the results for this
            particular trial.
        inherent_structure : bool
        """

        results_directory = Path(root) / Path("results")
        if inherent_structure:
            files = glob2.glob(str(results_directory / f"*_psi_config_IS.txt"))
        else:
            files = glob2.glob(str(results_directory / f"*_psi_config.txt"))
        res = read_files_via_numpy(files)  # Returns a list

        # Construct dictionaries out of each result while also keeping
        # track of the max key value. There's probably a way to do this
        # using list comprehension, but this is easier for now, and the
        # total number of results will never more more than ~10k, so it
        # should be ok.
        dict_res = []
        max_key = 0
        for res_arr in res:
            d = {
                int(key): int(value) for key, value
                in zip(res_arr[:, 0], res_arr[:, 1])
            }
            for key, value in d.items():
                if key > max_key and value > 0:
                    max_key = key
            dict_res.append(d)

        # Evaluate the statistics
        stats = np.zeros(shape=(len(dict_res), max_key + 1))
        for ii, d in enumerate(dict_res):
            for key, value in d.items():
                if value > 0:
                    stats[ii, key] = value

        # Get the location of the directory to save the finalized results
        if inherent_structure:
            final_path = Path(root) / Path("final/psi_config_IS.txt")
        else:
            final_path = Path(root) / Path("final/psi_config.txt")

        grid = [2**ii for ii in range(stats.shape[1])]
        N = stats.shape[0]
        psi = stats.mean(axis=0).squeeze()
        norm = psi.sum()

        # Compute psi as a "probability distribution"; normalize the sd and
        # standard errors as well.
        psi = psi / norm
        psi_sd = stats.std(axis=0).squeeze() / norm
        psi_sterr = psi_sd / np.sqrt(N - 1)

        to_save = np.array([grid, psi, psi_sd, psi_sterr])
        np.savetxt(final_path, to_save.T)  # Save as a numpy binary

        print(f"\tPsi Config (IS={inherent_structure}) done")

    @staticmethod
    def aging_config(root):
        """Evalutes all saved pi/aging config results. Evaluates both the
        standard trajectory and inherent structure results.

        Notes
        -----
        The loaded grids have the form:
        COLUMN 0: grid
        COLUMN 1: config index
        COLUMN 2: integer representation of the config
        COLUMN 3: integer representation of the inherent structure config

        Parameters
        ----------
        root : str
            The base directory containing all of the results for this
            particular trial.
        """

        # Load the data. Note it's really important to sort the file lists each
        # time so that the line up properly in the eq == pi1 == pi2 element
        # wise boolean comparison.
        results_directory = Path(root) / Path("results")
        files1 = glob2.glob(str(results_directory / f"*_pi1_config.txt"))
        files1.sort()
        pi1 = np.array(read_files_via_numpy(files1))
        files2 = glob2.glob(str(results_directory / f"*_pi2_config.txt"))
        files2.sort()
        pi2 = np.array(read_files_via_numpy(files2))

        # Need the pi grid too. This is contained in pi1 but we'll load it
        # separately here just to be sure.
        pi_grid = np.loadtxt(Path(root) / Path("grids/pi1.txt"))

        # Find where the config index, the standard config integer rep, and
        # the inherent structure rep are equal. This object, eq, is of shape
        # e.g. (1000, 91, 4) == (trials, recorded grid, columns in notes
        # above).
        eq = pi1 == pi2

        # Get the statistics
        mu = eq.mean(axis=0)
        sd = eq.std(axis=0)
        N = eq.shape[0]
        stderr = sd / np.sqrt(N - 1)

        # Parse through the columns and save each.
        results_index = np.array([pi_grid, mu[:, 1], sd[:, 1], stderr[:, 1]])
        results_stand = np.array([pi_grid, mu[:, 2], sd[:, 2], stderr[:, 2]])
        results_IS = np.array([pi_grid, mu[:, 3], sd[:, 3], stderr[:, 3]])

        final_path = Path(root) / Path("final/aging_config_index.txt")
        np.savetxt(final_path, results_index.T)  # Save as a numpy binary

        final_path = Path(root) / Path("final/aging_config.txt")
        np.savetxt(final_path, results_stand.T)

        final_path = Path(root) / Path("final/aging_config_IS.txt")
        np.savetxt(final_path, results_IS.T)

        print(f"\tAging Config done")

    @staticmethod
    def get_general_filename(
        base_fname, inherent_structure, energetic_threshold,
        extra_text=None, extension=".txt"
    ):
        """[summary]

        [description]

        Parameters
        ----------
        base_fname : str
            The base filename, such as "_psi_basin".
        inherent_structure : bool
        energetic_threshold : bool, optional
            If None, will ignore appending any string corresponding to if this
            is the energetic threshold or not.
        extra_text : str, optional
        extension : str, optional
            The etension to append at the end of the final string (the default
            is ".txt").

        Returns
        -------
        str
            The full file name.
        """

        f = base_fname

        if energetic_threshold is not None:
            if energetic_threshold:
                f += "_E"
            else:
                f += "_S"

        if inherent_structure:
            f += "_IS"

        if extra_text is not None:
            assert isinstance(extra_text, str)
            f += extra_text

        f += extension

        return f

    @staticmethod
    def psi_basin(root, inherent_structure, energetic_threshold):
        """Evaluates all psi basin results.

        Parameters
        ----------
        root : str
            The base directory containing all of the results for this
            particular trial.
        inherent_structure : bool
        energetic_threshold : bool
        """

        results_directory = Path(root) / Path("results")

        fname = Evaluator.get_general_filename(
            "*_psi_basin", inherent_structure, energetic_threshold,
            extra_text=None, extension=".txt"
        )
        files = glob2.glob(str(results_directory / Path(fname)))
        res = read_files_via_numpy(files)  # Returns a list

        # Construct dictionaries out of each result while also keeping
        # track of the max key value. There's probably a way to do this
        # using list comprehension, but this is easier for now, and the
        # total number of results will never more more than ~10k, so it
        # should be ok.
        dict_res = []
        dict_res_unique_per_basin = []
        max_key = 0
        for res_arr in res:

            # psi basin proper
            d = {
                int(key): int(value) for key, value
                in zip(res_arr[:, 0], res_arr[:, 1])
            }
            for key, value in d.items():
                if key > max_key and value > 0:
                    max_key = key
            dict_res.append(d)

            # Number of unique configs per basin
            d = {
                int(key): int(value) for key, value
                in zip(res_arr[:, 0], res_arr[:, 2])
            }
            for key, value in d.items():
                if key > max_key and value > 0:
                    max_key = key
            dict_res_unique_per_basin.append(d)

        # psi basin proper
        stats = np.zeros(shape=(len(dict_res), max_key + 1))
        for ii, d in enumerate(dict_res):
            for key, value in d.items():
                if value > 0:
                    stats[ii, key] = value

        # We only save results if the results if there is more than one entry
        # in psi. If there isn't, the most likely explanation is that the
        # barrier was either extremely high or ifinitely large/undefined
        # such as in the EREM case with beta = 2.
        if stats.shape[1] > 1:

            # Get the statistics
            grid = [2**ii for ii in range(stats.shape[1])]
            mu = stats.mean(axis=0)
            norm = mu.sum()
            mu = mu / norm
            N = stats.shape[0]
            sd = stats.std(axis=0) / norm
            stderr = sd / np.sqrt(N - 1)

            fname = Evaluator.get_general_filename(
                "psi_basin", inherent_structure, energetic_threshold,
                extra_text=None, extension=".txt"
            )
            final_path = Path(root) / Path(f"final/{fname}")
            stats_final = np.array([grid, mu.squeeze(), sd, stderr])
            np.savetxt(final_path, stats_final.T)

        # unique configs per basin
        stats_u = np.zeros(shape=(len(dict_res_unique_per_basin), max_key + 1))
        for ii, d in enumerate(dict_res_unique_per_basin):
            for key, value in d.items():
                if value > 0:
                    stats_u[ii, key] = value

        if stats_u.shape[1] > 1:

            # Get the statistics
            grid = [2**ii for ii in range(stats_u.shape[1])]
            mu = stats_u.mean(axis=0)
            N = stats_u.shape[0]
            sd = stats_u.std(axis=0)
            stderr = sd / np.sqrt(N - 1)

            fname = Evaluator.get_general_filename(
                "psi_basin", inherent_structure, energetic_threshold,
                extra_text="_unique_configs_per", extension=".txt"
            )
            final_path = Path(root) / Path(f"final/{fname}")
            stats_u_final = np.array([grid, mu.squeeze(), sd, stderr])
            np.savetxt(final_path, stats_u_final.T)

        print(
            f"\tPsi Basin (IS={inherent_structure}, "
            f"E={energetic_threshold}) done"
        )

    @staticmethod
    def _aging_basin_helper(arr, inherent_structure, energetic_threshold):
        """Helper function for the aging_basin calculation.

        Computes the meaningful statistics of the provided array. Handles the
        cases of inherent structure or which energetic threshold to use.

        Parameters
        ----------
        arr : list
            A list containing the two arrays [arr1, arr2].
        inherent_structure : bool
        energetic_threshold : bool

        Returns
        -------
        mu, sd, stderr : np.array
            Numpy arrays for the mean, standard deviation and standard errors
            of the results.
        """

        if not inherent_structure and energetic_threshold:
            slice1 = 1
            slice2 = 2
        elif inherent_structure and energetic_threshold:
            slice1 = 3
            slice2 = 4
        elif not inherent_structure and not energetic_threshold:
            slice1 = 5
            slice2 = 6
        elif inherent_structure and not energetic_threshold:
            slice1 = 7
            slice2 = 8
        else:
            raise NotImplementedError

        mu = []
        sd = []
        npts = []

        for gridpoint in range(arr[0].shape[1]):

            # This indexes whether or not the tracer is in a basin at pi1.
            # We only consider these points
            current_arr_pi1 = arr[0][:, gridpoint, slice2]

            in_basin = np.where(current_arr_pi1 == 1)[0]
            to_consider_pi1 = arr[0][in_basin, gridpoint, slice1]
            to_consider_pi2 = arr[1][in_basin, gridpoint, slice1]

            # We only care about when those basin indexes are equal
            eq = to_consider_pi1 == to_consider_pi2
            npts.append(eq.shape[0])
            mu.append(np.mean(eq))
            sd.append(np.std(eq))

        mu = np.array(mu)
        sd = np.array(sd)
        npts = np.array(npts)  # Not a scalar!

        return mu, sd, sd / np.sqrt(npts - 1)

    @staticmethod
    def aging_basin(root):
        """Evalutes all saved pi/aging config results.

        Notes
        -----
        Column 0: grid point
        ENERGY:
        Column 1: basin index
        Column 2: in basin or not (1 if in basin, 0 if not in basin)
        Column 3: IS basin index
        Column 4: IS in basin or not
        ENTROPY:
        Column 5: basin index
        Column 6: in basin or not (1 if in basin, 0 if not in basin)
        Column 7: IS basin index
        Column 8: IS in basin or not
        """

        results_directory = Path(root) / Path("results")
        files1 = glob2.glob(str(results_directory / f"*_pi1_basin.txt"))
        files1.sort()
        pi1 = np.array(read_files_via_numpy(files1))
        files2 = glob2.glob(str(results_directory / f"*_pi2_basin.txt"))
        files2.sort()
        pi2 = np.array(read_files_via_numpy(files2))

        # Need the pi grid too. This is contained in pi1 but we'll load it
        # separately here just to be sure.
        pi_grid = np.loadtxt(Path(root) / Path("grids/pi1.txt"))

        combinations = [
            (True, True), (False, False), (True, False), (False, True)
        ]
        for (b1, b2) in combinations:
            mu, sd, stderr = Evaluator._aging_basin_helper([pi1, pi2], b1, b2)
            fname = Evaluator.get_general_filename(
                "aging_basin", b1, b2, extension=".txt"
            )
            final_path = Path(root) / Path(f"final/{fname}")
            to_save = np.array([pi_grid, mu, sd, stderr])
            np.savetxt(final_path, to_save.T)

            print(f"\tAging Basin (IS={b1}, E={b2}) done")

    @staticmethod
    def _process_ridge_energy(arr):
        """Processes a single array's worth of results for the ridge energy
        calculation. We only keep data which is representative of a simulation
        where the tracer has actually gone from below a threshold to above it.
        There are quite a few simulations where this never occurs. Default
        values for mu and sd in these cases are 0.0, and the default values for
        the maxima and minima are -1e15 and 1e15, respectively, obviously
        posing a problem if these data points were included in the averaging.
        This method removes those data points.

        Parameters
        ----------
        arr : np.array
            A numpy array of shape (ntracers, 6).
        """

        keep = np.where(
            (arr[:, 1] != 0.0) & (arr[:, 2] > 0.0) & (arr[:, 0] == 1)
        )
        arr = arr[keep]

        # Collect averaged quantities:
        # Starting with the weighted mean. We treat the ground truth value here
        # as the mean ridge energy weighted by the number of data points
        # contributing to that average.

        # See the comments by the except ValueError below.
        if arr.shape[0] < 2:
            return None

        try:
            weighted = arr[:, 1] * arr[:, 5]
            mu = np.mean(weighted)
            sd = np.std(weighted)
            sderr = sd / np.sqrt(weighted.shape[0] - 1)
            total_max = np.max(arr[:, 3])  # Max of max
            total_min = np.min(arr[:, 4])  # Min of min
            return np.array([mu, sd, sderr, total_max, total_min])

        # Catches the value error "zero-size array to reduction operation
        # maximum which has no identity". All this is saying is that there are
        # no data points in the array after keeping only the sensible ones.
        # In this case, we won't save any results, since this indicates no
        # valid data.
        except ValueError:
            return None

    @staticmethod
    def ridge_energy(root, inherent_structure, energetic_threshold):
        """Evaluates the ridge energy information.

        res has the shape of e.g. (ntracers, 2, 6). The second dimension
        indexes whether or not the energy right before jumping above a
        threshold is the same as right after it (0, 1, respectively). The
        third dimension indexes the recorded quantity:
            0: Either 0 or 1. Indexes whether or not the tracer ever exited the
               first basin. If this is 0, that data point should be ignored.
            1: Average ridge energy value. If this value is 0, this also
               indicates that this data point should be ignored.
            2: Variance of the ridge energy values.
            3: The maximum ridge energy value.
            4: The minimum ridge energy value.
            5: The total number of points that have contributed to the average.
        """

        results_directory = Path(root) / Path("results")
        fname = Evaluator.get_general_filename(
            "*_ridge", inherent_structure, energetic_threshold,
            extra_text=None, extension=".txt"
        )
        files = glob2.glob(str(results_directory / Path(fname)))
        files.sort()
        res = np.array(read_files_via_numpy(files))  # See docstring

        res_same = Evaluator._process_ridge_energy(res[:, 0, :].squeeze())
        if res_same is not None:
            fname = Evaluator.get_general_filename(
                "ridge_energy", inherent_structure, energetic_threshold,
                extra_text="_same", extension=".txt"
            )
            final_path = Path(root) / Path(f"final/{fname}")
            np.savetxt(final_path, res_same[:, None].T)

        res_diff = Evaluator._process_ridge_energy(res[:, 1, :].squeeze())
        if res_diff is not None:
            fname = Evaluator.get_general_filename(
                "ridge_energy", inherent_structure, energetic_threshold,
                extra_text="_diff", extension=".txt"
            )
            final_path = Path(root) / Path(f"final/{fname}")
            np.savetxt(final_path, res_diff[:, None].T)

        print(
            f"\tRidge Energy (IS={inherent_structure}, "
            f"E={energetic_threshold}) done"
        )

    def __init__(self, specified_directory=None):
        cache = u.get_cache()
        _all_trial_dirs = u.listdir_fp(cache)
        all_trial_dirs = [d for d in _all_trial_dirs if os.path.isdir(d)]
        if specified_directory is not None:
            all_trial_dirs = [
                d for d in all_trial_dirs if specified_directory in d
            ]
        print(f"Evaluating {len(all_trial_dirs)} directories...")
        for full_dir_path in all_trial_dirs:
            print(f"Evaluating {Path(full_dir_path).stem}")
            Evaluator.energy(full_dir_path)
            Evaluator.psi_config(full_dir_path, inherent_structure=False)
            Evaluator.psi_config(full_dir_path, inherent_structure=True)
            Evaluator.aging_config(full_dir_path)
            Evaluator.psi_basin(full_dir_path, True, True)
            Evaluator.psi_basin(full_dir_path, False, False)
            Evaluator.psi_basin(full_dir_path, False, True)
            Evaluator.psi_basin(full_dir_path, True, False)
            Evaluator.aging_basin(full_dir_path)
            Evaluator.ridge_energy(full_dir_path, True, True)
            Evaluator.ridge_energy(full_dir_path, False, False)
            Evaluator.ridge_energy(full_dir_path, False, True)
            Evaluator.ridge_energy(full_dir_path, True, False)
        print("All done")
