# hdspin

Sandbox for the Exponential and Gaussian Random Energy Models. If you use this code, please consider citing our [work](https://doi.org/10.1103/PhysRevE.106.024603).

## Installation instructions

The hdspin repository requires no external libraries whatsoever (except via Git Submodules), (mostly) everything is self-contained. Downloading and building the code should be simple via CMake. Note, it is required you have MPI available on your system in order to build hdspin. You can do this via your system's package managers (such as Homebrew or apt). hdspin is tested with openmpi. See [here](https://github.com/mpi4py/setup-mpi/blob/master/setup-mpi.sh) for how hdspin's CI system installs MPI (you can emulate this).

```bash
git clone git@github.com:matthewcarbone/hdspin.git
cd hdspin
git submodule init
git submodule update
cmake -S . -B build -DPRECISON=256 -DBUILD_TESTS=ON -DSMOKE=ON
cd build
make
```

There are three options for the user to set:
* `-DPRECISON=<INT>` is the maximum number of spins you can use during the simulation. Should be a power of 2 (as recommended by the Arbitrary Precision library). Default is `256`. Note that this is the _maximum_ value you can use for `N_spins` in the simulation.
* `-DBUILD_TESTS={ON, OFF}` is a boolean flag for telling CMake whether or not to compile the testing suite. Default is `OFF`.
* `-DSMOKE={ON, OFF}` controls whether or not to use the smoke testing or not. Smoke tests basically run tests using slightly less statistics, and are generally faster. Default is `ON`.

## Running instructions

### Running the hdspin executable

After running `make` in the previous steps, an executable `build/hdspin` will be created. Running `hdspin` is simple. `hdspin` uses a command line parser called [CLI11](https://github.com/CLIUtils/CLI11). Use `hdspin -h` to see a list of options. A `config.json` is always saved to the working directory with all of the command line inputs.

Four parameters are absolutely required:
* `log10_N_timesteps=<INT>`: the log10 number of timesteps to run 
* `N_spins=<INT>`: the number of spins to use in the simulation. Must be `<=PRECISON`.
* `beta=<FLOAT>`: inverse temperature (`beta_critical` is set automatically based on the `landscape`).
* `landscape={"EREM", "GREM"}`: the type of simulation to run (either exponential or Gaussian REM).


### Post-processing

Once simulations are computed, run the post-processing script from the same working directory to finalize the results into the `final` directory

```bash
python3 /path/to/postprocess.py
```

Post-processing creates averages and spreads of all observable quantities, such as the energy.

# License

The `hdspin` code is released under a 3-clause BSD license. Hosted codes are contained locally as per the permissive terms of the associated licenses. This includes nlohmann's [Json](https://github.com/nlohmann/json) header, as well as [Catch2](https://github.com/catchorg/Catch2), [CLI11](https://github.com/CLIUtils/CLI11) and the [Arbitrary Precision](https://www.codeproject.com/Articles/5319814/Arbitrary-Precision-Easy-to-use-Cplusplus-Library) library.


# Funding acknowledgement

This software is based upon work supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Department of Energy Computational Science Graduate Fellowship under Award Number DE-FG02-97ER25308.
