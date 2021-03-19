#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import numpy as np
import os
import pickle


class PlottingManager:

    def __init__(
        self, capsize=2, capthick=0.3, elw=0.3, marker='s', ms=1.0,
        lw=1.0
    ):
        self.plot_kwargs = {
            'linewidth': lw,
            'marker': marker,
            'ms': ms,
            'capthick': capthick,
            'capsize': capsize,
            'elinewidth': elw
        }

    def plot_energy(
        self, ax, directory, cache=os.environ['HDSPIN_CACHE_DIR'],
        fname='final/energy.txt', standard_error=False, color=None,
        inherent_structure=False, label=None
    ):
        """Plotting for the energy.

        Parameters
        ----------
        ax
            The matplotlib axis to use.
        directory : str
            The location of the directory to use.
        cache : str
            The cache location, defaults to HDSPIN_CACHE_DIR environment
            variable.
        fname : str
            The name of the energy file.
        standard_error : bool
            Whether or not to use the standard error vs standard deviation for
            plotting errorbars.
        """

        arr = np.loadtxt(os.path.join(cache, directory, fname))
        grid = arr[:, 0]

        div = 1.0
        if standard_error:
            div = np.sqrt(arr.shape[0] - 1)

        if not inherent_structure:
            e = arr[:, 1]
            e_sd = arr[:, 2]

        else:
            e = arr[:, 3]
            e_sd = arr[:, 4]

        ax.errorbar(
            grid, e, yerr=e_sd / div, color=color, label=label,
            **self.plot_kwargs
        )

    def plot_psi_config(
        self, ax, directory, cache=os.environ['HDSPIN_CACHE_DIR'],
        fname_base='final/psi_config', standard_error=False, color=None,
        label=None, inherent_structure=False
    ):

        if inherent_structure:
            fname = fname_base + "_IS.txt"
        else:
            fname = fname_base + ".txt"
        arr = np.loadtxt(os.path.join(cache, directory, fname))
        grid = [2**ii for ii in range(arr.shape[1])]
        e = arr.mean(axis=0)
        e_sd = arr.std(axis=0)
        div = 1.0
        if standard_error:
            div = np.sqrt(arr.shape[0] - 1)
        norm = e.sum()
        ax.errorbar(
            grid, e / norm, yerr=e_sd / norm / div, color=color, label=label,
            **self.plot_kwargs
        )

    def plot_aging_config(
        self, ax, directory, cache=os.environ['HDSPIN_CACHE_DIR'],
        fname='final/aging_config_sd.pkl', standard_error=False, color=None,
        label=None, ctype='index', grid_loc='grids/pi1.txt'
    ):
        """Choose ctype from index, standard, and inherent_structure."""

        arr = pickle.load(open(os.path.join(cache, directory, fname), 'rb'))
        pi_grid = np.loadtxt(os.path.join(cache, directory, grid_loc))

        if ctype == 'index':
            slice = 1
        elif ctype == 'standard':
            slice = 2
        elif ctype == 'inherent_structure':
            slice = 3
        else:
            raise NotImplementedError

        pi_val = arr[:, :, slice]
        mu = pi_val.mean(axis=0)
        sd = pi_val.std(axis=0)

        div = 1.0
        if standard_error:
            div = np.sqrt(arr.shape[0] - 1)

        ax.errorbar(
            pi_grid, mu, yerr=sd / div, color=color, label=label,
            **self.plot_kwargs
        )

    def plot_aging_basin(
        self, ax, directory, cache=os.environ['HDSPIN_CACHE_DIR'],
        fname='final/aging_basin_sd.pkl', standard_error=False, color=None,
        label=None, ctype='standard', grid_loc='grids/pi1.txt', threshold='E',
        div_x_axis_by=1.0
    ):
        """Choose ctype from index, standard, and inherent_structure."""

        arr = pickle.load(open(os.path.join(cache, directory, fname), 'rb'))
        pi_grid = np.loadtxt(os.path.join(cache, directory, grid_loc))

        if ctype == 'standard' and threshold == 'E':
            slice1 = 1  # Basin index
            slice2 = 2  # In basin (1) or not (0)
        elif ctype == 'inherent_structure' and threshold == 'E':
            slice1 = 3
            slice2 = 4
        elif ctype == 'standard' and threshold == 'S':
            slice1 = 5
            slice2 = 6
        elif ctype == 'inherent_structure' and threshold == 'S':
            slice1 = 7
            slice2 = 8
        else:
            raise NotImplementedError

        mu = []
        sd = []
        npts = []
        # print(arr[0].shape)
        for gridpoint in range(arr[0].shape[1]):

            # This indexes whether or not the tracer is in a basin at pi1.
            # We only consider these points
            current_arr_pi1 = arr[0][:, gridpoint, slice2]
            # current_arr_pi2 = arr[1][:, gridpoint, slice2]
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
        npts = np.array(npts)

        div = np.ones(shape=(npts.shape[0]))
        if standard_error:
            div = np.sqrt(npts - 1)

        ax.errorbar(
            pi_grid / div_x_axis_by, mu, yerr=sd / div, color=color,
            label=label, **self.plot_kwargs
        )

    def plot_psi_basin(
        self, ax, directory, cache=os.environ['HDSPIN_CACHE_DIR'],
        fname_base='final/psi_basin', standard_error=False, color=None,
        label=None, inherent_structure=False, threshold='E', style='-',
        unique_configs=False
    ):

        if inherent_structure:
            fname = fname_base + f"_{threshold}_IS"
        else:
            fname = fname_base + f"_{threshold}"
        if unique_configs:
            fname = fname + "_u"
        fname = fname + ".txt"
        arr = np.loadtxt(os.path.join(cache, directory, fname))
        grid = [2**ii for ii in range(arr.shape[1])]
        e = arr.mean(axis=0)
        e_sd = arr.std(axis=0)
        div = 1.0
        if standard_error:
            div = np.sqrt(arr.shape[0] - 1)
        norm = e.sum()
        ax.errorbar(
            grid, e / norm, yerr=e_sd / norm / div, color=color, label=label,
            linestyle=style, **self.plot_kwargs
        )
