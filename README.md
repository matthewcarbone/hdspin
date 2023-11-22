<div align=center>

# hdspin

**Lightning-fast simulator for the Exponential and Gaussian Random Energy Models** <br>

_If you use this code, please consider citing our [work](https://doi.org/10.1103/PhysRevE.106.024603)_ <br> 

---

[![image](https://github.com/matthewcarbone/hdspin/actions/workflows/build.yaml/badge.svg)](https://github.com/matthewcarbone/hdspin/actions/workflows/build.yaml)
<!-- [![image](https://github.com/matthewcarbone/hdspin/actions/workflows/tests.yml/badge.svg)](https://github.com/matthewcarbone/hdspin/actions/workflows/tests.yml) -->

</div>

# 🚀 Features

⚡ Fast spin-flips using [decimal representation and bit-flip operations](https://github.com/matthewcarbone/hdspin/blob/master/inc/spin.h), driven by the [Arbitrary Precision](https://www.codeproject.com/Articles/5319814/Arbitrary-Precision-Easy-to-use-Cplusplus-Library) (AP) library, no need to store an expensive `std::vector<int>` for every system anymore.

⚡ System energies are stored using a [least recently used cache](https://www.geeksforgeeks.org/lru-cache-implementation/), leading to orders of magnitude more memory saving and together wth the AP library, arbitrarily large systems.

⚡ [MPI load-balancer](https://github.com/matthewcarbone/hdspin/blob/master/src/main_utils.cpp), allowing for massively parallel simulations on high-performance computing systems. The rare time-consuming job no longer holds up other simulations.

# 📕 Installation instructions

hdspin requires [MPI](http://www.mpi-forum.org), and is tested using [MPICH](https://www.mpich.org) and [Open MPI](https://www.open-mpi.org) on both Ubuntu latest and MacOS. Other than that, every external dependency is self-contained explicitly under the terms of their licences. Installing hdspin should be straightforward using CMake:

```bash
git clone git@github.com:matthewcarbone/hdspin.git
cd hdspin
cmake -S . -B build -DPRECISON=256 -DBUILD_TESTS=ON -DSMOKE=ON
cd build
make
```

At compile time, there are three options for the user to set:
* `-DPRECISON=<INT>` is the maximum number of spins you can use during the simulation. Should be a power of 2 (as recommended by the Arbitrary Precision library). Default is `256`. Note that this is the _maximum_ value you can use for `N_spins` in the simulation.
* `-DBUILD_TESTS={ON, OFF}` is a boolean flag for telling CMake whether or not to compile the testing suite. Default is `OFF`.
* `-DSMOKE={ON, OFF}` controls whether or not to use the smoke testing or not. Smoke tests basically run tests using slightly less statistics, and are generally faster. Default is `ON`.

# 📗 Running instructions

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


# 📘 Post-processing

Once simulations are computed, run the post-processing script from the same working directory to finalize the results into the `final` directory

```bash
python3 /path/to/postprocess.py
```

Post-processing creates averages and spreads of all observable quantities, such as the energy.

# 🔨 License

The hdspin software is released under a 3-clause BSD license.


# 🙏 Software acknowledgement

Much of this work would have been significantly harder if not for some incredible products of the open source community! 

The hdspin software hosts codes locally as per the permissive terms of the associated licenses (check out this explicitly clear [blog post](https://levelofindirection.com/blog/unit-testing-in-cpp-and-objective-c-just-got-ridiculously-easier-still.html) regarding the Catch library for more details, if you're interested). These excellent packages are:

- Niels Lohmann's [Json](https://github.com/nlohmann/json) header, for handling config serialization to disk
- Alexander Ponomarev's [LRU cache](https://github.com/lamerman/cpp-lru-cache) header, for handling the... LRU cache
- [Catch2](https://github.com/catchorg/Catch2) (the header-only version), for unit testing
- [CLI11](https://github.com/CLIUtils/CLI11) for all command line argument needs
- The [Arbitrary Precision](https://www.codeproject.com/Articles/5319814/Arbitrary-Precision-Easy-to-use-Cplusplus-Library) library, for all arbitrary precision integer uses

# 💲 Funding acknowledgement

This software is based upon work supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Department of Energy Computational Science Graduate Fellowship under Award Number DE-FG02-97ER25308.
