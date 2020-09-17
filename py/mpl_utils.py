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
