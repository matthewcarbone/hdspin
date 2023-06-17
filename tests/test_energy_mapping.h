#include "energy_mapping.h"
#include "utils.h"
#include "utils_testing_suite.h"

bool _test_energy_mapping_sampling_EREM_given_beta_critical(const double beta_critical)
{
    parameters::SimulationParameters sp;
    sp.landscape = "EREM";
    sp.beta_critical = beta_critical;
    sp.use_manual_seed = true;
    sp.seed = 1234;
    EnergyMapping emap = EnergyMapping(sp);
    const double mean = -1.0 / sp.beta_critical;
    const double variance = 1.0 / sp.beta_critical / sp.beta_critical;
    const int N = 1000000;
    const double eps = 0.1;
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


bool _test_energy_mapping_sampling_REM_given_N_spins(const int N_spins)
{
    parameters::SimulationParameters sp;
    sp.landscape = "REM";
    sp.N_spins = N_spins;
    sp.use_manual_seed = true;
    sp.seed = 4567;
    EnergyMapping emap = EnergyMapping(sp);
    const double mean = 0.0;
    const double variance = N_spins;
    const int N = 10000000;
    const double eps = 0.1;
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
