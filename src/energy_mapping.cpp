#include <cstring>
#include <random>

#include "energy_mapping.h"
#include "utils.h"
#include "lru.h"
#include "utils.h"
#include "ArbitraryPrecision/ap/ap.hpp"


double EnergyMapping::sample_energy() const
{
    if (params.landscape == "EREM")
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

void EnergyMapping::get_config_energies_array_(const ap_uint<PRECISON> *neighbors, double *neighboring_energies, const unsigned int bitLength) const
{
    for (int ii=0; ii<bitLength; ii++)
    {
        neighboring_energies[ii] = get_config_energy(neighbors[ii]);
    }
}

void EnergyMapping::_initialize_distributions()
{
    if (params.use_manual_seed == true)
    {
        generator.seed(params.seed);
    }
    else
    {
        const unsigned int seed = std::random_device{}();  
        generator.seed(seed);
    }
    
    // Initialize the distributions themselves
    // If the distribution type is not found throws a runtime_error
    if (params.landscape == "EREM")
    {
        const double p = params.beta_critical;
        exponential_distribution.param(
            std::exponential_distribution<double>::param_type(p)
        );
    }
    else if (params.landscape == "GREM")
    {
        const double p = sqrt(params.N_spins);
        normal_distribution.param(
            std::normal_distribution<double>::param_type(0.0, p)
        );
    }
    else
    {
        const std::string err = "Invalid landscape " + params.landscape;
        throw std::runtime_error(err);
    }
}


EnergyMapping::EnergyMapping(const parameters::SimulationParameters params) : params(params)
{
    _initialize_distributions();

    // Allocate the energy storage mediums
    if (params.memory == -1){
        const long long n_configs = pow(2, params.N_spins);
        energy_map.set_capacity(n_configs);
    }
    else if (params.memory > 0){energy_map.set_capacity(params.memory);}
    else
    {
        throw std::runtime_error("Invalid choice for memory; must be either -1 or >0");
    }
};


unsigned int _min_element(const double *arr, const unsigned int N)
{
    unsigned int min_el = 0;
    for (int ii=1; ii<N; ii++)
    {
        if (arr[ii] < arr[min_el]){min_el = ii;}
    }
    return min_el;
}

ap_uint<PRECISON> EnergyMapping::get_inherent_structure(const ap_uint<PRECISON> state) const
{
    unsigned int min_el;
    double tmp_energy;
    ap_uint<PRECISON> tmp_state = state;

    ap_uint<PRECISON>* tmp_neighbors = 0;
    tmp_neighbors = new ap_uint<PRECISON> [params.N_spins];

    double* tmp_neighbor_energies = 0;
    tmp_neighbor_energies = new double [params.N_spins];

    while (true)
    {
        state::get_neighbors_(tmp_neighbors, tmp_state, params.N_spins);
        get_config_energies_array_(tmp_neighbors, tmp_neighbor_energies, params.N_spins);
        min_el = _min_element(tmp_neighbor_energies, params.N_spins);
        tmp_energy = get_config_energy(tmp_state);

        if (tmp_neighbor_energies[min_el] < tmp_energy)
        {
            // flip to new energy
            tmp_state = tmp_neighbors[min_el];
        }
        else{break;}
    }

    delete[] tmp_neighbors;
    delete[] tmp_neighbor_energies;

    return tmp_state;
}
