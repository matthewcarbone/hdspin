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

    # This is for the mean median
    weights = arr[:, :, 2]
    dat = arr[:, :, 0]
    mu = np.ma.average(dat, axis=0, weights=weights)
    var = np.ma.average((mu - dat)**2, axis=0, weights=weights)

    # This is for the mean mean
    weights = arr[:, :, 2]
    dat = arr[:, :, 1]
    mu2 = np.ma.average(dat, axis=0, weights=weights)
    var2 = np.ma.average((mu - dat)**2, axis=0, weights=weights)

    final = np.array([
        grid,
        mu,
        np.sqrt(var),
        np.sqrt(var) / np.sqrt(weights.sum(axis=0)),
        mu2,
        np.sqrt(var2),
        np.sqrt(var2) / np.sqrt(weights.sum(axis=0))
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


# def aging_config(all_filenames, substring, save_path):
#     """Process all the aging config results."""

#     files = [f for f in all_filenames if substring in f]
#     if len(files) == 0:
#         return
#     arr = np.array(read_files_via_numpy(files))
#     eq = arr[:, :, 0] == arr[:, :, 1]
#     mu = eq.mean(axis=0)
#     sd = eq.std(axis=0)
#     stderr = sem(eq, axis=0)
#     # grid = np.loadtxt("grids/pi1.txt")
#     grid = np.array([2**ii for ii in range(arr.shape[1])])
#     final = np.array([grid, mu, sd, stderr]).T
#     np.savetxt(save_path, final, fmt='%.08e')


if __name__ == '__main__':
    Path(FINAL_DIRECTORY).mkdir(exist_ok=True, parents=False)
    fnames = get_all_results_filenames()

    # Handle all standard one-point observables
    obs1(fnames, "_energy.txt", Path(FINAL_DIRECTORY) / "energy.txt")
    obs1(fnames, "_energy_IS.txt", Path(FINAL_DIRECTORY) / "energy_IS.txt")
    ridge(fnames, "_ridge_E.txt", Path(FINAL_DIRECTORY) / "ridge_E.txt")
    ridge(fnames, "_ridge_S.txt", Path(FINAL_DIRECTORY) / "ridge_S.txt")
    obs1(fnames, "_acceptance_rate.txt", Path(FINAL_DIRECTORY) / "acceptance_rate.txt")
    obs1(fnames, "_inherent_structure_timings.txt", Path(FINAL_DIRECTORY) / "inherent_structure_timings.txt")
    obs1(fnames, "_walltime_per_waitingtime.txt", Path(FINAL_DIRECTORY) / "walltime_per_waitingtime.txt")

    # All two point observables
    psi_config(fnames, "_psi_config.txt", Path(FINAL_DIRECTORY) / "psi_config.txt")
    psi_basin(fnames, "_psi_basin_E.txt", Path(FINAL_DIRECTORY) / "psi_basin_E.txt")
    psi_basin(fnames, "_psi_basin_S.txt", Path(FINAL_DIRECTORY) / "psi_basin_S.txt")

    # ...
    cache_size(fnames, "_cache_size.txt", Path(FINAL_DIRECTORY) / "cache_size.txt")
