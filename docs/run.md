# Running instructions

After running `make` in the previous steps, an executable `build/hdspin` will be created. Running hdspin is simple. hdspin uses a command line parser called [CLI11](https://github.com/CLIUtils/CLI11). Use `hdspin -h` to see a list of options. A `config.json` is always saved to the working directory with all of the command line inputs. All outputs are saved in the working directory as well.

Four parameters are absolutely required:
* `log10_N_timesteps=<INT>`: the log10 number of timesteps to run 
* `N_spins=<INT>`: the number of spins to use in the simulation. Must be `<=PRECISON`.
* `beta=<FLOAT>`: inverse temperature (`beta_critical` is set automatically based on the `landscape`).
* `landscape={"EREM", "GREM"}`: the type of simulation to run (either exponential or Gaussian REM).

One example of a job might be

```bash
mpiexec -n 5 path/to/hdspin -N 20 -l EREM -b 2.5 -t 6 -n 100 --seed=123
```

which will run the exponential random energy model with 20 spins with inverse temperature `beta=2.5`, for `1e6` timesteps and 100 tracers (with seed 123 for reproducibility). The job will be split amongst 4 compute tasks with a single controller task (for 5 total).


# Post-processing

Results are post-processed automatically to `results.json` at the end of the simulation.
