# hdspin
REM/EREM sandbox

## Installation instructions

The `hdspin` repository requires no external libraries whatsoever, everything is self-contained. Building the code should be simple via CMake.

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

After running `make` in the previous steps, an executable `build/hdspin` will be created. Running `hdspin` is simple. `hdspin` takes a single optional argument: the path to the `config.json` file to be used. If no argument is provided, then it is assumed that `config.json` is in the current working directory. Note that all results will be saved to the current working directory.

```bash
mpiexec -n <N_TASKS> /path/to/hdspin [optional config path]
```

The `config.json` file absolutely requires 4 parameters (the rest are set to default values):
* `log10_N_timesteps <INT>`: the log10 number of timesteps to run 
* `N_spins <INT>`: the number of spins to use in the simulation. Must be `<=PRECISON`.
* `beta <FLOAT>`: inverse temperature (`beta_critical` is set automatically based on the `landscape`).
* `landscape {"EREM", "GREM"}`: the type of simulation to run (either exponential or Gaussian REM).

```json
{
    "log10_N_timesteps": 5,
    "N_spins": 20,
    "beta": 2.8,
    "landscape" : "EREM",
}
```

Once simulations are computed, run the postprocessing script from the same working directory to finalize the results into the `final` directory

```bash
python3 /path/to/postprocess.py
```

## License

The `hdspin` code is released under a 3-clause BSD license. Hosted codes are contained locally as per the permissive terms of the associated licenses. This includes nlohmann's [Json](https://github.com/nlohmann/json) header, as well as [Catch2](https://github.com/catchorg/Catch2) and the [Arbitrary Precision](https://www.codeproject.com/Articles/5319814/Arbitrary-Precision-Easy-to-use-Cplusplus-Library) library.


## Funding acknowledgement

This software is based upon work supported by the
* U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Department of Energy Computational Science Graduate Fellowship under Award Number DE-FG02-97ER25308;
*  U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences under Contract No. DE-SC0012704; 
* U.S. Department of Energy, Office of Science, Office of Workforce Development for Teachers and Scientists (WDTS) under the Visiting Faculty Program (VFP).

The Software resulted from work developed under a U.S. Government
Contract No. DE-SC0012704 and are subject to the following terms: the
U.S. Government is granted for itself and others acting on its behalf a
paid-up, nonexclusive, irrevocable worldwide license in this computer
software and data to reproduce, prepare derivative works, and perform
publicly and display publicly.

THE SOFTWARE IS SUPPLIED \"AS IS\" WITHOUT WARRANTY OF ANY KIND. THE
UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND THEIR
EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE, TITLE OR NON-INFRINGEMENT, (2) DO NOT ASSUME
ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF THE
SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4) DO NOT WARRANT
THAT THE SOFTWARE WILL FUNCTION UNINTERRUPTED, THAT IT IS ERROR-FREE OR
THAT ANY ERRORS WILL BE CORRECTED.

IN NO EVENT SHALL THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
ENERGY, OR THEIR EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF ANY KIND OR
NATURE RESULTING FROM EXERCISE OF THIS LICENSE AGREEMENT OR THE USE OF
THE SOFTWARE.





















