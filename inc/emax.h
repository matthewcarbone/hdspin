#include "spin.h"
#include "utils.h"

#include <map>

struct EMaxt2Data
{
    double max_energy = -1e16;
    double tmin;
    double tmax;
};

class EMaxt2
{
private:
    const utils::SimulationParameters params;
    const SpinSystem* spin_system_ptr;
    size_t length;
    mutable std::map<size_t, EMaxt2Data> trackers;
    mutable std::map<size_t, double> max_energies;

public:
    EMaxt2(const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double simulation_clock) const;
    json as_json() const;
};
