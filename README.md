# hdspin
REM/EREM sandbox

## Installation instructions

The `hdspin` repository requires no external libraries whatsoever, everything is self-contained. Building the code should be simple via CMake.

```bash
git clone git@github.com:matthewcarbone/hdspin.git
cd hdspin
# Checkout to appropriate branch,
# probably master but could be mc/remove-all-vectors right now 
git checkout mc/remove-all-vectors
git submodule init
git submodule update
cmake -S . -B build -DPRECISON=256 -DBUILD_TESTS=ON -DSMOKE=ON 
cd build
make
```

There are three options for the user to set:
* `-DPRECISON=<INT>` is the maximum number of spins you can use during the simulation. Should be a power of 2 (as recommended by the Arbitrary Precision library). Default is `256`.
* `-DBUILD_TESTS={ON, OFF}` is a boolean flag for telling CMake whether or not to compile the testing suite. Default is `OFF`.
* `-DSMOKE={ON, OFF}` controls whether or not to use the smoke testing or not. Smoke tests basically run tests using slightly less statistics, and are generally faster. Default is `ON`.

## License

The `hdspin` code is released under a 3-clause BSD license. Hosted codes are contained locally as per the permissive terms of the associated licenses. This includes nlohmann's [Json](https://github.com/nlohmann/json) header, as well as [Catch2](https://github.com/catchorg/Catch2) and the [Arbitrary Precision](https://www.codeproject.com/Articles/5319814/Arbitrary-Precision-Easy-to-use-Cplusplus-Library) library.
























