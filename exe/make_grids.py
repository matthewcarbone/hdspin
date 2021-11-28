import sys

import numpy as np


def make_grids(dw, log10_timesteps, energy_gridpoints, pi_gridpoints):
    """Writes the grids to disk in the standard spot."""

    nMC = int(10**log10_timesteps)

    energy_grid = np.unique(np.logspace(
        0, np.log10(nMC), energy_gridpoints, dtype=int, endpoint=True
    ))

    # Adds 0 to the energy grid so we can record the value of the energy
    # at 0 itself.
    energy_grid = np.array([0] + list(energy_grid))

    tw_max = nMC // (dw + 1.0)

    # Define the first grid.
    pi_g1 = np.unique(np.logspace(
        0, np.log10(tw_max), pi_gridpoints, dtype=int, endpoint=True
    ))

    # The second grid is directly related to the first via
    # tw -> tw + tw * dw
    pi_g2 = (pi_g1 * (dw + 1.0)).astype(int)

    np.savetxt("grids/energy.txt", energy_grid, fmt="%i")
    np.savetxt("grids/pi1.txt", pi_g1, fmt="%i")
    np.savetxt("grids/pi2.txt", pi_g2, fmt="%i")


if __name__ == '__main__':
    make_grids(
        float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
        int(sys.argv[4])
    )
