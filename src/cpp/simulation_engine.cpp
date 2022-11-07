#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/stl.h"
#include <iostream>
#include <vector>
#include <algorithm>  // sorting

namespace py = pybind11;


/**
 * @brief Container for the results of the simulation.
 * @details Contains multiple useful objects that are created during the
 * running of any of the simulations. These include a vector of n_record times
 * and n_record int-vectors containing the states at those times.
 * 
 */
struct SimulationResult
{
    int test_return = 0;
    int test_return2 = 0;
    std::vector<std::vector<int> > vector_configurations;
};


SimulationResult run_simulation(const int N, const int M)
{
    SimulationResult simulation_result;
    simulation_result.test_return = N;
    simulation_result.test_return2 = M;
    return simulation_result;
}


PYBIND11_MODULE(simulation_engine, m) {
    py::class_<SimulationResult>(m, "SimulationResult")
        .def_readwrite("test_return", &SimulationResult::test_return)
        .def_readwrite("test_return2", &SimulationResult::test_return2)
        .def_readwrite("vector_configurations", &SimulationResult::vector_configurations);
    m.def("run_simulation", &run_simulation, "Executes the spin simulation");
}
