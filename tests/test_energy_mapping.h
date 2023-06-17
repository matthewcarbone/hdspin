#include "energy_mapping.h"
#include "utils.h"

bool _test_energy_mapping_sampling_EREM_given_beta_critical(const double beta_critical)
{
    parameters::SimulationParameters sp;
    sp.landscape = "EREM";
    sp.beta_critical = beta_critical;
    sp.use_manual_seed = true;
    sp.seed = 1234;
    EnergyMapping emap = EnergyMapping(sp);
    const double mean = -1.0 / sp.beta_critical;
    double v = 0.0;
    const int N = 1000000;
    const double eps = 0.01;
    for (int ii=0; ii<N; ii++)
    {
        v += emap.sample_energy();
    }
    v /= double(N);
    // std::cout << v << " " << mean << std::endl;
    if ((v < mean + eps) && (v > mean - eps)){return true;}
    return false;
}
