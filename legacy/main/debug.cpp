#include <fstream>      // std::ofstream
#include <chrono>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <sstream>
#include <omp.h>

#include "Utils/structures.h"
#include "Spin/gillespie.h"
#include "Spin/standard.h"
#include "Engine/sim.h"


void print_int_vector(const std::vector<int> p)
{ 
    for (int i = 0; i < p.size(); i++) 
    {
       std::cout << p[i] << " "; 
    }
    std::cout << std::endl;
}


int main(int argc, char *argv[])
{
    const int log_N_timesteps = 6;
    const int N_spins = 8;
    const double beta = 1.333;
    const double beta_critical = 1.0;
    const int landscape = 0;  // 0 for EREM, 1 for REM
    // const int dynamics = 1;  // 0 for standard, 1 for gillespie
    const RuntimeParameters params = get_runtime_params(log_N_timesteps,
        N_spins, beta, beta_critical, landscape);

    GillespieSpinSystem sys(params);
    std::cout << "s1" << std::endl;
    std::vector<int> c = sys.get_spin_config();
    print_int_vector(c);
    sys.flip_spin_(2);
    c = sys.get_spin_config();
    print_int_vector(c);

    FileNames fnames;
    StandardSimulation sim(fnames, params);
    sim.execute();
}
