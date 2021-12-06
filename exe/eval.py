#!/usr/bin/env python3


import os
from pathlib import Path
import sys
import warnings

import glob2
import numpy as np
from scipy.stats import sem


RESULTS_DIRECTORY = Path("data")


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


def get_all_results_filenames():
    """Gets a list of all of the filenames in the provided path

    Parameters
    ----------
    path : posix Path, optional

    Returns
    -------
    list
    """

    fnames = list(Path("data").iterdir())
    return [str(f) for f in fnames]


def energy(all_filenames, substring, save_path):
    """Processes files that have returned energy information."""

    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return
    arr = np.array(read_files_via_numpy(files))
    grid = np.loadtxt("grids/energy.txt")
    mu = arr.mean(axis=0)
    sd = arr.std(axis=0)
    stderr = sem(arr, axis=0)
    final = np.array([grid, mu, sd, stderr]).T
    np.savetxt(save_path, final, fmt='%.08e')


def ridges(all_filenames, substring, save_path):
    """Process the ridge energies by concatenating all of the results."""

    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        files = read_files_via_numpy(files)
    files = [np.atleast_2d(np.array(f)) for f in files if len(f) > 0]
    arr = np.concatenate(files, axis=0)
    np.savetxt(save_path, arr, fmt='%.08e')


def psi_config(all_filenames, substring, save_path):
    """Process all of the config information."""

    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return
    files = read_files_via_numpy(files)
    arr = np.array([np.array(f) for f in files if len(f) > 0])
    grid = np.array([2**ii for ii in range(arr.shape[1])])
    psi = arr.mean(axis=0)
    psi = psi / psi.sum()
    sd = arr.std(axis=0)
    stderr = sem(arr, axis=0)
    final = np.array([grid, psi, sd, stderr]).T
    where = np.where(final[:, 1] != 0)[0]
    final = final[where, :]
    np.savetxt(save_path, final, fmt='%.08e')


def aging_config(all_filenames, substring, save_path):
    """Process all the aging config results."""

    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return
    arr = np.array(read_files_via_numpy(files))
    eq = arr[:, :, 0] == arr[:, :, 1]
    mu = eq.mean(axis=0)
    sd = eq.std(axis=0)
    stderr = sem(eq, axis=0)
    grid = np.loadtxt("grids/pi1.txt")
    final = np.array([grid, mu, sd, stderr]).T
    np.savetxt(save_path, final, fmt='%.08e')


def psi_basin(all_filenames, substring, save_path):
    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return
    arr = np.array(read_files_via_numpy(files))
    arr = arr[:, :, 0].squeeze()  # ignore unique configs/basin
    grid = np.array([2**ii for ii in range(arr.shape[1])])
    mu = arr.mean(axis=0)
    sd = arr.std(axis=0)
    stderr = sem(arr, axis=0)
    final = np.array([grid, mu, sd, stderr]).T
    where = np.where(final[:, 1] != 0)[0]
    final = final[where, :]
    np.savetxt(save_path, final, fmt='%.08e')


def aging_basin(all_filenames, substring, save_path):
    """"""

    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        arr = np.array(read_files_via_numpy(files))

    mu = []
    sd = []
    stderr = []

    for gridpoint in range(arr.shape[1]):
        subarr = arr[:, gridpoint, :]

        # Indexes are vec idx 1, in basin 1, vec idx 2, in basin 2
        in_basin_1 = np.where(subarr[:, 1] == 1)[0]
        eq = subarr[in_basin_1, 0] == subarr[in_basin_1, 2]

        mu.append(np.mean(eq))
        sd.append(np.std(eq))
        stderr.append(sem(eq))

    grid = np.loadtxt("grids/pi1.txt")
    final = np.array([grid, mu, sd, stderr]).T
    np.savetxt(save_path, final, fmt='%.08e')


if __name__ == '__main__':
    Path("final").mkdir(exist_ok=True, parents=False)
    fnames = get_all_results_filenames()
    energy(fnames, "_energy.txt", "final/energy.txt")
    energy(fnames, "_energy_IS.txt", "final/energy_IS.txt")
    ridges(fnames, "ridge_E.txt", "final/ridge_E.txt")
    ridges(fnames, "ridge_S.txt", "final/ridge_S.txt")
    psi_config(fnames, "_psi_config.txt", "final/psi_config.txt")
    psi_config(fnames, "_psi_config_IS.txt", "final/psi_config_IS.txt")
    aging_config(fnames, "_aging_config.txt", "final/aging_config.txt")
    aging_config(fnames, "_aging_config_IS.txt", "final/aging_config_IS.txt")
    psi_basin(fnames, "_psi_basin_E.txt", "final/psi_basin_E.txt")
    psi_basin(fnames, "_psi_basin_E_IS.txt", "final/psi_basin_E_IS.txt")
    psi_basin(fnames, "_psi_basin_S.txt", "final/psi_basin_S.txt")
    psi_basin(fnames, "_psi_basin_S_IS.txt", "final/psi_basin_S_IS.txt")
    aging_basin(fnames, "_aging_basin_E.txt", "final/aging_basin_E.txt")
    aging_basin(fnames, "_aging_basin_E_IS.txt", "final/aging_basin_E_IS.txt")
    aging_basin(fnames, "_aging_basin_S.txt", "final/aging_basin_S.txt")
    aging_basin(fnames, "_aging_basin_S_IS.txt", "final/aging_basin_S_IS.txt")
