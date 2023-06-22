#ifndef TEST_SPIN_H
#define TEST_SPIN_H

#include "spin.h"

bool test_basic_lru_cache_with_spin()
{
    parameters::SimulationParameters p;
    p.log10_N_timesteps = 4;
    p.N_timesteps = ipow(10, int(p.log10_N_timesteps));
    p.N_spins = 4;
    p.landscape = "EREM";
    p.beta = 2.4;
    p.beta_critical = 1.0;
    p.dynamics = "standard";
    p.memory = 15;
    p.n_tracers_per_MPI_rank = 1;

    EnergyMapping emap(p);
    SpinSystem spin(p, emap);
    std::vector<double> vec;

    for (int ii=0; ii<16; ii++)
    {
        spin.set_state(ii);
        // std::cout << spin.binary_state() << " " << spin.energy() << std::endl;
        vec.push_back(spin.energy());
    }
    for (int ii=15; ii>=0; ii--)
    {
        spin.set_state(ii);
        // std::cout << spin.binary_state() << " " << spin.energy() << std::endl;
        vec.push_back(spin.energy());
    }

    int ii1, ii2;
    for (int ii=0; ii<16; ii++)
    {
        ii1 = ii;
        ii2 = 32 - ii - 1;

        // Only the 0th state should be different
        // std::cout << ii1 << " " << ii2 << " " << vec[ii1] << " " << vec[ii2] << std::endl;
        if ((ii1 == 0) && (vec[ii1] == vec[ii2])){return false;}
        else if ((ii1 != 0) && (vec[ii1] != vec[ii2])){return false;}
    }

    return true;
}

#endif
