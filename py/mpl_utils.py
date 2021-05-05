#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import numpy as np
from pathlib import Path
from scipy.integrate import quad

from py.utils import listdir_fp


def hxw(x, w=0.5):
    """Equation 3 in the Cammarota 2018 paper."""

    def integrand(u):
        return 1.0 / (1.0 + u) / u**x

    return quad(integrand, w, np.inf)[0] * np.sin(np.pi * x) / np.pi


class AnalyticResults:
    """Manages the analytic results for the REM and EREM models.

    Parameters
    ----------
    model : {"EREM", "REM"}, optional
        The model type. (The default is "EREM").
    """

    def __init__(self, model="EREM"):
        self.model = model

    def energy(
        self, ax, g, e, beta, color='gray', linestyle='--', label=None,
        anchor=-1, zorder=3, offset=0.5
    ):
        """Plots the correct slope of the energy, on a logarithmic x-scale.

        Parameters
        ----------
        ax
            The matplotlib axis to plot on.
        g, e : array_like
            The time and energy grids.
        beta : float
            The inverse temperature.
        color, linestyle, label : str
        anchor : int
            The point on the grids to anchor the analytic line to.
        zorder : int
        offset : float
            The y-offset of the line to be plotted.
        """

        if self.model == "EREM":
            slope = -1.0 / beta
            b = e[anchor] - slope * np.log(g[anchor])
            ax.plot(
                g, b + slope * np.log(g) + offset, color=color, linestyle='--',
                label=label, zorder=zorder
            )
        else:
            raise NotImplementedError

    def psi(
        self, ax, g, beta, color='gray', linestyle='--', label=None,
        anchor=-1, zorder=3, offset=0.5
    ):

        if self.model == "EREM":
            ax.plot(
                g, offset*g**(-1.0 * 1.0 / beta), color=color, linestyle='--',
                label=label, zorder=zorder
            )
        else:
            raise NotImplementedError

    def aging(
        self, ax, g0, gf, beta, w=0.5, color='black', label=None,
        linestyle='--', zorder=3
    ):
        """[summary]

        [description]
        """

        h = hxw(1.0 / beta, w=w)
        ax.hlines(
            h, g0, gf, color=color, label=label, linestyle=linestyle,
            zorder=zorder
        )


class ResultsManager:
    """Manages results and plotting.

    Parameters
    ----------
    directory : str
        The location of the results. Directories should look something like
        e.g., abs/path/to/standard_erem_8_24_1333_20210322221923
    """

    def __init__(
        self,
        directory,
        capsize=2,
        capthick=0.3,
        elw=0.3,
        marker='s',
        ms=1.0,
        lw=1.0
    ):

        # Set some parameters for uniform plotting.
        self.plot_kwargs = {
            'linewidth': lw,
            'marker': marker,
            'ms': ms,
            'capthick': capthick,
            'capsize': capsize,
            'elinewidth': elw
        }

        # Load in the results
        self.results = {
            Path(f).stem: np.loadtxt(f, delimiter=" ")
            for f in listdir_fp(directory)
        }

    @staticmethod
    def get_key(
        result, inherent_structure=None, energetic_threshold=None,
        unique_configs_per=None, diff=None
    ):
        """Gets the key for self.results easily.

        Parameters
        ----------
        result : {'energy', 'energy_eg_traj'}
        inherent_structure, energetic_threshold : bool

        Returns
        -------
        str
            The key for self.results.
        """

        if result == "energy":
            if inherent_structure:
                return "energy_IS"
            return "energy"
        elif result == "energy_eg_traj":
            if inherent_structure:
                return "energy_IS_eg_traj"
            return "energy_eg_traj"
        elif result == "psi_basin":
            k = "psi_basin"
            if energetic_threshold:
                k += "_E"
            else:
                k += "_S"
            if inherent_structure:
                k += "_IS"
            if unique_configs_per:
                k += "_unique_configs_per"
            return k
        elif result == "aging_basin":
            k = "aging_basin"
            if energetic_threshold:
                k += "_E"
            else:
                k += "_S"
            if inherent_structure:
                k += "_IS"
            return k
        elif result == "ridge_energy":
            k = "ridge_energy"
            if energetic_threshold:
                k += "_E"
            else:
                k += "_S"
            if inherent_structure == 1:
                k += "_IS"
            elif inherent_structure == 2:
                k += "_proxy_IS"
            if diff:
                k += "_diff"
            else:
                k += "_same"
            return k

        return None

    def _plot_standard(
        self, ax, arr, standard_error, color, label, linestyle
    ):
        """Plots an array (arr) on the provided axis.

        Parameters
        ----------
        arr : np.array
            A numpy array of the specific shape (N, 4). The first column is
            the x axis, the second, the y axis, the third is the standard
            deviation and fourth is standard error.
        standard_error : bool
            If True, plots the 4th column. If False, plots the 3rd.
        """

        ax.errorbar(
            arr[:, 0], arr[:, 1],
            yerr=arr[:, 3] if standard_error else arr[:, 2],
            color=color, label=label, linestyle=linestyle, **self.plot_kwargs
        )

    def energy(
        self, ax=None, plot=False, standard_error=True, color='black',
        label=None, inherent_structure=False, single_trajectory=0,
        linestyle="-"
    ):
        """Plotting for the energy.

        Parameters
        ----------
        ax
            The matplotlib axis to use.
        standard_error : bool
            Whether or not to use the standard error vs standard deviation for
            plotting errorbars. (The default is True).
        color : str
            The color of the plot. (The default is 'black').
        label : str
            The legend label for the plot. (The default is None).
        inherent_structure : bool
            If True, loads the inherent structure trajectories instead of the
            standard ones. (The default is False).
        single_trajectory : int
            If 0, will plot the average energy with errorbars. If > 0, will
            plot that many trajectories with random colors. Note that this
            ignores the standard_error, color and label arguments.

        Returns
        -------
        np.array
        """

        if single_trajectory == 0:
            arr = self.results[ResultsManager.get_key(
                "energy", inherent_structure=inherent_structure
            )]

            if plot:
                self._plot_standard(
                    ax, arr, standard_error, color, label,
                    linestyle
                )
        else:
            arr = self.results[ResultsManager.get_key(
                "energy_eg_traj", inherent_structure=inherent_structure
            )]

            if plot:
                for ii in range(single_trajectory):
                    ax.errorbar(
                        arr[:, 0], arr[:, ii + 1], yerr=None,
                        linestype=linestyle, **self.plot_kwargs
                    )

        return arr

    def psi_config(
        self, ax=None, plot=False, standard_error=True, color='black',
        label=None, inherent_structure=False
    ):
        """See the `energy` method for details."""

        arr = self.results["psi_config_IS"] if inherent_structure \
            else self.results["psi_config"]

        if plot:
            self._plot_standard(ax, arr, standard_error, color, label)

        return arr

    def psi_basin(
        self, ax=None, plot=False, standard_error=True, color='black',
        label=None, inherent_structure=False, energetic_threshold=True,
        unique_configs_per=False, linestyle="-"
    ):
        """See the `energy` method for details."""

        key = ResultsManager.get_key(
            "psi_basin", inherent_structure, energetic_threshold,
            unique_configs_per
        )

        try:
            arr = self.results[key]
        except KeyError:
            return None

        if plot:
            self._plot_standard(
                ax, arr, standard_error, color, label, linestyle
            )

        return arr

    def aging_config(
        self, ax=None, plot=False, standard_error=True, color='black',
        label=None, res_type='standard', linestyle="-"
    ):
        """See the `energy` method for details.

        Parameters
        ----------
        res_type : {'standard', 'IS', 'index'}
        """

        assert res_type in ['standard', 'IS', 'index']

        key = "aging_config"
        if res_type == "IS":
            key += "_IS"
        elif res_type == "index":
            key += "_index"
        arr = self.results[key]

        if plot:
            self._plot_standard(
                ax, arr, standard_error, color, label, linestyle
            )

        return arr

    def aging_basin(
        self, ax=None, plot=False, standard_error=True, color='black',
        label=None, inherent_structure=False, energetic_threshold=True,
        linestyle="-"
    ):
        """See the `energy` method for details.

        Parameters
        ----------
        res_type : {'standard', 'IS', 'index'}
        """

        key = ResultsManager.get_key(
            "aging_basin", inherent_structure=inherent_structure,
            energetic_threshold=energetic_threshold
        )
        arr = self.results[key]

        if plot:
            self._plot_standard(
                ax, arr, standard_error, color, label, linestyle
            )

        return arr

    def ridge_energy(
        self, inherent_structure=0, energetic_threshold=True, diff=False
    ):
        """[summary]

        Plotting must be done manually.
        [mu, sd, sderr, total_max, total_min, ...hist..., ...bins...]
        """

        key = ResultsManager.get_key(
            "ridge_energy", inherent_structure=inherent_structure,
            energetic_threshold=energetic_threshold, diff=diff
        )
        try:
            res = self.results[key]
        except KeyError:
            return None

        return res
