# Installation instructions

hdspin requires [CMake](https://cmake.org) to build and [MPI](http://www.mpi-forum.org) to run. It is tested using [MPICH](https://www.mpich.org) and [Open MPI](https://www.open-mpi.org) on both Ubuntu and MacOS latest. Otherwise, every external dependency is self-contained explicitly under the terms of their licences. Installing hdspin should be straightforward using CMake:

## Install CMake and MPI

There are a variety of ways to install CMake and MPI. For example, on Mac, I prefer to use Homebrew:

```bash
brew install cmake
brew install open-mpi
```

Or on Linux, something like this (depending on the distribution you're using and its particular package manager)

```bash
sudo apt-get update
sudo apt-get install cmake
sudo apt install build-essential
sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
```

> [!WARNING]
> hdspin is not tested on Windows, please use at your own risk!


## Install hdspin

Once CMake and MPI are installed on your system, the rest should be straighforward. First, download the hdspin source code (I prefer `git clone` but you can also download a tarball). `cd` into `hdspin` and use CMake to generate all required makefiles. Then run make (with optional `j8` to use 8 parallel jobs) Something like this should work just fine:

```bash
git clone git@github.com:matthewcarbone/hdspin.git
cd hdspin
cmake -S . -B build -DPRECISON=256 -DBUILD_TESTS=ON -DSMOKE=ON
cd build
make -j8
```

Note that at compile time, there are three options for the user to set:
* `-DPRECISON=<INT>` is the maximum number of spins you can use during the simulation. Should be a power of 2 (as recommended by the Arbitrary Precision library). Default is `256`. Note that this is the _maximum_ value you can use for `N_spins` in the simulation.
* `-DBUILD_TESTS={ON, OFF}` is a boolean flag for telling CMake whether or not to compile the testing suite. Default is `OFF`.
* `-DSMOKE={ON, OFF}` controls whether or not to use the smoke testing or not. Smoke tests basically run tests using slightly less statistics, and are generally faster. Default is `ON`.
