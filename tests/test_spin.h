#ifndef TEST_SPIN_H
#define TEST_SPIN_H

#include "spin.h"
#include "utils.h"

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


// bool test_inherent_structure()
// {
//     parameters::SimulationParameters p;
//     p.log10_N_timesteps = 4;
//     p.N_timesteps = ipow(10, int(p.log10_N_timesteps));
//     p.N_spins = 4;
//     p.landscape = "EREM";
//     p.beta = 2.4;
//     p.beta_critical = 1.0;
//     p.dynamics = "standard";
//     p.memory = -1;
//     p.n_tracers_per_MPI_rank = 1;
//     p.use_manual_seed = true;
//     p.seed = 123;

//     EnergyMapping emap(p);
//     SpinSystem sys(p, emap);

//     for (ap_uint<PRECISON> state=0; state<16; state++)
//     {
//         sys.set_state(state);
//         std::cout << state << " "  << sys.binary_state() << " " << sys.energy() << " " << sys.get_inherent_structure() << std::endl;
//     }
//     return true;

// }

bool test_inherent_structure_min_is_min()
{
    const unsigned int N = 12;
    parameters::SimulationParameters p;
    p.log10_N_timesteps = 4;
    p.N_timesteps = ipow(10, int(p.log10_N_timesteps));
    p.N_spins = N;
    p.landscape = "EREM";
    p.beta = 2.4;
    p.beta_critical = 1.0;
    p.dynamics = "standard";
    p.memory = -1;
    p.n_tracers_per_MPI_rank = 1;
    p.use_manual_seed = true;
    p.seed = 123;

    EnergyMapping emap(p);
    SpinSystem sys(p, emap);
    EnergyMapping* emap_ptr = &emap;

    ap_uint<PRECISON>* neighbors = 0;
    neighbors = new ap_uint<PRECISON> [p.N_spins];

    double* neighboring_energies = 0;
    neighboring_energies = new double [p.N_spins];

    const unsigned int m = pow(2, N);
    for (ap_uint<PRECISON> s=0; s<m; s++)
    {
        sys.set_state(s);
        // ap_uint<PRECISON> inherent_structure = sys._help_get_inherent_structure();

        ap_uint<PRECISON> inherent_structure = emap_ptr->get_inherent_structure(s);

        // std::cout << s << " " << inherent_structure << std::endl;

        // If the s is the inherent structure, run the local
        // assertion
        if (s == inherent_structure)
        {
            state::get_neighbors_(neighbors, s, p.N_spins);
            emap.get_config_energies_array_(neighbors, neighboring_energies, p.N_spins);
            for (int jj=0; jj<p.N_spins; jj++)
            {
                if (neighboring_energies[jj] < sys.energy()){return false;}
            }
            // std::cout << s << " "  << sys.binary_state() << " " << sys.energy() << " " << sys.get_inherent_structure() << std::endl;
        }
        
    }

    delete[] neighbors;
    delete[] neighboring_energies;

    return true;
}


#endif
