<div align=center>
  
# What does hdspin calculate?

</div>

**We outline the various observables that hdspin calculates.**

> [!IMPORTANT]
> Statistics are not saved at every timestep. As it is possible to have a _huge_ number of timesteps, possibly on the order of many billions, this is not tractable from a storage standpoint. Grids are generally linear in log-space, and observables are saved only at those timesteps. Most observables have a variety of statistics saved, including the mean, standard deviation, standard error, and median (at every timestep on the grid).

All averaged statistics can be found in the `results.json` file saved at the end of the simulation. You can load the results via the built-in `json` library in Python via something like this:
```python3
import json
with open("results.json", "r") as f:
  results = json.load(f)
```
Similarly, diagnostics are saved in `diagnostics.json` and the original config is saved in `config.json`. These files can be loaded simiarly. You can access the available observables via `results.keys()`, and for a given observable, such as the energy, you can access the available statistical information via `results["energy"].keys()`.

# Energy

Consider the energy of tracer $i$ at simulation clock time $t$ is $E_{i}(t).$ The average energy over $N$ tracers is simply given by

$$ E(t) = \frac{1}{N} \sum_{i=1}^N E_i(t).$$

Results for the energy can be accessed via `results["energy"]`. Similarly, the grid on which the energy is recorded can be accessed via `results["grids"]["energy"]`. The standard deviation, and standard error, are similarly calculated and accessed,

$$ \sigma_E(t) = \sqrt{\frac{1}{N} \sum_{i=1}^N (E(t) - E_i(t))^2}.$$

All in summary, one can access these variables via:

```python
key = "energy"
grid = results["grids"][key]
energy = results[key]["mean"]
energy_standard_error = results[key]["standard_error"]
```

# E-max

TK

# Ridge energy

# Autocorrelation functions ($\psi$)

TK

# Aging functions ($\Pi$)

TK

# Walltime-per-waiting time

TK

# Acceptance rate

The cumulative acceptance rate is logged on the same grid as the energy. It represents the Metropolis acceptance rate up until that point in time, and can be accessed via `key = "acceptance_rate"`. Defining the cumulative acceptance rate at time $t$ for tracer $i$ as $A_i(t),$ it is given by the ratio of the number of accepted steps to the number of total steps. Note that in Gillespie simulations, $A_i(t) = 1.$ The statistical moments of the cumulative acceptance rate are calculated in the same fashion as the energy.

# Cache size

TK



