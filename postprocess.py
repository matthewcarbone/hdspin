#!/usr/bin/env python3

from pathlib import Path

import numpy as np
from scipy.stats import sem


RESULTS_DIRECTORY = Path("data")
FINAL_DIRECTORY = Path("final")


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

    fnames = list(Path(RESULTS_DIRECTORY).iterdir())
    return [str(f) for f in fnames]


def obs1(all_filenames, substring, save_path):
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


def ridge(all_filenames, substring, save_path):
    """Processes files that have returned energy information."""

    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return
    arr = np.array(read_files_via_numpy(files))
    grid = np.loadtxt("grids/energy.txt")
    mu = np.ma.average(arr[:, :, 0], axis=0, weights=arr[:, :, 1])
    var = np.ma.average((mu - arr[:, :, 0])**2, axis=0, weights=arr[:, :, 1])
    final = np.array([
        grid, mu, np.sqrt(var),
        np.sqrt(var) / arr[:, :, 1].sum(axis=0)
    ]).T
    np.savetxt(save_path, final, fmt='%.08e')


def cache_size(all_filenames, substring, save_path):

    files = [f for f in all_filenames if substring in f]
    if len(files) == 0:
        return
    arr = np.array(read_files_via_numpy(files))
    n = arr[0, 0]
    arr = arr[:, 1:] / n  # Percentage of cache filled
    grid = np.loadtxt("grids/energy.txt")
    mu = arr.mean(axis=0)
    sd = arr.std(axis=0)
    stderr = sem(arr, axis=0)
    final = np.array([grid, mu, sd, stderr]).T
    np.savetxt(save_path, final, fmt='%.08e')


if __name__ == '__main__':
    Path(FINAL_DIRECTORY).mkdir(exist_ok=True, parents=False)
    fnames = get_all_results_filenames()

    # Handle all _standard_ one-point observables
    obs1(fnames, "_energy.txt", Path(FINAL_DIRECTORY) / "energy.txt")
    obs1(fnames, "_energy_IS.txt", Path(FINAL_DIRECTORY) / "energy_IS.txt")
    ridge(fnames, "_ridge_E.txt", Path(FINAL_DIRECTORY) / "ridge_E.txt")
    ridge(fnames, "_ridge_S.txt", Path(FINAL_DIRECTORY) / "ridge_S.txt")
    obs1(fnames, "_acceptance_rate.txt", Path(FINAL_DIRECTORY) / "acceptance_rate.txt")
    obs1(fnames, "_inherent_structure_timings.txt", Path(FINAL_DIRECTORY) / "inherent_structure_timings.txt")
    obs1(fnames, "_walltime_per_waitingtime.txt", Path(FINAL_DIRECTORY) / "walltime_per_waitingtime.txt")

    # ...
    cache_size(fnames, "_cache_size.txt", Path(FINAL_DIRECTORY) / "cache_size.txt")
