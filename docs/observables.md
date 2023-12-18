# Observables

We outline the various observables that hdspin calculates.

- Statistics are not saved at every timestep. As it is possible to have many millions of timesteps, this is not tractable from a storage standpoint. Grids are generally linear in log-space, and observables are saved only at those timesteps.
- Most observables have a variety of statistics saved, including the mean, standard deviation, standard error, and median (at every timestep on the grid).  
- All averaged statistics can be found in the `results.json` file saved at the end of the simulation.

## Average energy

Consider the energy of tracer $i$ at simulation clock time $t$ is $E_{i}(t).$ The average energy over $N$ tracers is simply given by

$$ E(t) = \frac{1}{N} \sum_{i=1}^N E_i(t).$$

Results for the energy can be accessed via `results["energy"]`.
