#ifndef TEST_ENERGY_MAPPING_H
#define TEST_ENERGY_MAPPING_H

#include "energy_mapping.h"
#include "utils.h"
#include "utils_testing_suite.h"

namespace test_energy_mapping
{

bool test_energy_mapping_sampling_EREM_given_beta_critical(const double beta_critical)
{
    parameters::SimulationParameters sp;
    sp.landscape = "EREM";
    sp.beta_critical = beta_critical;
    sp.use_manual_seed = true;
    sp.seed = 1234;
    EnergyMapping emap = EnergyMapping(sp);
    const double mean = -1.0 / sp.beta_critical;
    const double variance = 1.0 / sp.beta_critical / sp.beta_critical;
    int N;
    double eps;
    if (SMOKE == 1)
    {
        N = 100000;
        eps = 1.0;
    }
    else
    {
        N = 1000000;
        eps = 0.1;    
    }
    std::vector<double> v;
    for (int ii=0; ii<N; ii++)
    {
        v.push_back(emap.sample_energy());
    }

    const double num_mean = _mean_vector(v);
    const double num_var = _variance_vector(v);

    // std::cout << num_mean << " " << mean << std::endl;
    // std::cout << num_var << " " << variance << std::endl; 

    if ((num_mean > mean + eps) || (num_mean < mean - eps)){return false;}    
    if ((num_var > variance + eps) || (num_var < variance - eps)){return false;}

    return true;
}


bool test_energy_mapping_sampling_REM_given_N_spins(const int N_spins)
{
    parameters::SimulationParameters sp;
    sp.landscape = "GREM";
    sp.N_spins = N_spins;
    sp.use_manual_seed = true;
    sp.seed = 4567;
    EnergyMapping emap = EnergyMapping(sp);
    const double mean = 0.0;
    const double variance = N_spins;

    int N;
    double eps;
    if (SMOKE == 1)
    {
        N = 1000000;
        eps = 1.0;
    }
    else
    {
        N = 10000000;
        eps = 0.1;    
    }
    
    std::vector<double> v;
    for (int ii=0; ii<N; ii++)
    {
        v.push_back(emap.sample_energy());
    }

    const double num_mean = _mean_vector(v);
    const double num_var = _variance_vector(v);

    // std::cout << num_mean << " " << mean << std::endl;
    // std::cout << num_var << " " << variance << std::endl;

    if ((num_mean > mean + eps) || (num_mean < mean - eps)){return false;}    
    if ((num_var > variance + eps) || (num_var < variance - eps)){return false;}

    return true;
}


bool test_small_cache(const int N_spins)
{
    parameters::SimulationParameters sp;
    sp.landscape = "GREM";
    sp.N_spins = N_spins;
    sp.use_manual_seed = true;
    sp.seed = 4567;
    sp.memory = 3;
    EnergyMapping emap = EnergyMapping(sp);
    
    const ap_uint<PRECISON> state_1 = 1234;
    const ap_uint<PRECISON> state_2 = 5678;
    const ap_uint<PRECISON> state_3 = 9123;
    const ap_uint<PRECISON> state_4 = 3456;

    const double e1 = emap.get_config_energy(state_1);
    // std::cout << e1 << " " << emap.get_config_energy(state_1) << std::endl;
    if (e1 != emap.get_config_energy(state_1)){return false;}

    const double e2 = emap.get_config_energy(state_2);
    // std::cout << e2 << " " << emap.get_config_energy(state_2) << std::endl;
    if (e2 != emap.get_config_energy(state_2)){return false;}

    const double e3 = emap.get_config_energy(state_3);
    // std::cout << e3 << " " << emap.get_config_energy(state_3) << std::endl;
    if (e3 != emap.get_config_energy(state_3)){return false;}

    const double e4 = emap.get_config_energy(state_4);
    // std::cout << e4 << " " << emap.get_config_energy(state_4) << std::endl;
    if (e4 != emap.get_config_energy(state_4)){return false;}

    bool same = (e1 == emap.get_config_energy(state_1));
    // std::cout << "should be different: " << e1 << " " << emap.get_config_energy(state_1) << " " << same << std::endl;
    if (same){return false;}

    same = (e4 == emap.get_config_energy(state_4));
    // std::cout << "should be the same: " << e4 << " " << emap.get_config_energy(state_4) << " " << same << std::endl;
    if (!same){return false;}

    // std::cout << " ----- " << std::endl;
    return true;
}


bool test_memory_minus_one(const int N_spins)
{
    parameters::SimulationParameters sp;
    sp.landscape = "GREM";
    sp.N_spins = N_spins;
    sp.use_manual_seed = true;
    sp.seed = 4567;
    sp.memory = -1;
    EnergyMapping emap = EnergyMapping(sp);

    // For every configuration, we call the get_config_energy()
    // function twice to make sure that nothing was popped.
    double e;
    for (int ii=0; ii<pow(2, N_spins); ii++)
    {
        e = emap.get_config_energy(ii);
        if (e != emap.get_config_energy(ii)){return false;}
    }
    return true;
}


bool test_massive_AP_LRU(const int N_spins)
{
    parameters::SimulationParameters sp;
    sp.landscape = "GREM";
    sp.N_spins = N_spins;
    sp.use_manual_seed = true;
    sp.seed = 45678;
    sp.memory = 1000;
    EnergyMapping emap = EnergyMapping(sp);

    ap_uint<PRECISON> large_number = arbitrary_precision_integer_pow(2, N_spins-1);

    double e;
    for (int ii=0; ii<10; ii++)
    {
        // Check that the mapping works each time
        // indicating that the lru cache is looking up the result
        // correctly
        e = emap.get_config_energy(large_number);
        if (e != emap.get_config_energy(large_number)){return false;}
        large_number -= 13;
    }
    return true;
}

}

#endif
