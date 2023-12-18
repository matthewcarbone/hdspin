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

# 📕 Documentation

Please see our documentation for instructions on how to [install](https://github.com/matthewcarbone/hdspin/blob/master/docs/install.md), [run](https://github.com/matthewcarbone/hdspin/blob/master/docs/run.md) and [understand the results](https://github.com/matthewcarbone/hdspin/blob/master/docs/observables.md) of hdspin!


# 🔨 License

The hdspin software is released under a [3-clause BSD license](https://github.com/matthewcarbone/hdspin/blob/master/LICENSE).


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
