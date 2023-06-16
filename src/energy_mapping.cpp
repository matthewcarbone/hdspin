#include <cstring>
#include <random>

#include "energy_mapping.h"
#include "structures.h"
#include "utils.h"
#include "lru.h"
#include "state_manipulation.h"
#include "ArbitraryPrecision/ap/ap.hpp"


double EnergyMapping::sample_energy() const
{
    if (rtp.landscape == "EREM")
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
    // Adjusted memory dynamics, only caching some
    // finite number of energies
    if (rtp.memory != 0)
    {
        // If our key exists in the LRU cache, simply return the value
        if (energy_map.key_exists(state))
        {
            return energy_map.get_fast(state);
        }

        // Otherwise, we sample a new value
        else
        {
            const double sampled = sample_energy();
            energy_map.put(state, sampled);
            return sampled;
        }
    }

    else{throw std::runtime_error(
        "Do not call get_config_energy in memoryless system");}
}
