#include <cstring>
#include <random>

#include "energy_mapping.h"
#include "utils.h"
#include "lru.h"
#include "utils.h"
#include "ArbitraryPrecision/ap/ap.hpp"


double EnergyMapping::sample_energy() const
{
    if (parameters.landscape == "EREM")
    {
        return -exponential_distribution(generator);
    }
    else
    {
        return normal_distribution(generator);
    }
}



double EnergyMapping::get_config_energy(const ap_uint<PRECISON> state) const
{

    const std::string state_string = std::string(state);

    // If our key exists in the LRU cache, simply return the value
    if (energy_map.key_exists(state_string))
    {
        return energy_map.get_fast(state_string);
    }

    // Otherwise, we sample a new value
    else
    {
        const double sampled = sample_energy();
        energy_map.put(state_string, sampled);
        return sampled;
    }
}


void EnergyMapping::_initialize_distributions()
{
    if (parameters.use_manual_seed == true)
    {
        generator.seed(parameters.seed);
    }
    else
    {
        const unsigned int seed = std::random_device{}();  
        generator.seed(seed);  
    }
    
    // Initialize the distributions themselves
    // If the distribution type is not found throws a runtime_error
    if (parameters.landscape == "EREM")
    {
        const double p = parameters.beta_critical;
        exponential_distribution.param(
            std::exponential_distribution<double>::param_type(p)
        );
    }
    else if (parameters.landscape == "REM")
    {
        const double p = sqrt(parameters.N_spins);
        normal_distribution.param(
            std::normal_distribution<double>::param_type(0.0, p)
        );
    }
    else
    {
        const std::string err = "Invalid landscape " + parameters.landscape;
        throw std::runtime_error(err);
    }
}


EnergyMapping::EnergyMapping(const parameters::SimulationParameters parameters) : parameters(parameters)
{
    _initialize_distributions();

    // Allocate the energy storage mediums
    if (parameters.memory == -1){
        const long long n_configs = pow(2, parameters.N_spins);
        energy_map.set_capacity(n_configs);
    }
    else if (parameters.memory > 0){energy_map.set_capacity(parameters.memory);}
    else
    {
        throw std::runtime_error("Invalid choice for memory; must be either -1 or >0");
    }
};
