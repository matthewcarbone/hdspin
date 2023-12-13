#ifndef OBS1_H
#define OBS1_H

#include <queue>
#include <vector>
#include <unordered_set>

#include "spin.h"
#include "utils.h"


class StreamingMedian
{
protected:

    // Standard  max priority queue
    std::priority_queue <double> max_heap;

    // Reverse (min) priority queue
    std::priority_queue <double, std::vector <double>, std::greater<double>> min_heap;

public:
    StreamingMedian();
    double median() const;
    void update(const double v);
};

class StreamingMean
{
protected:
    double value = 0.0;
    double counts = 0.0;

public:
    StreamingMean(){};
    double mean() const
    {
        if (counts > 0.0)
        {
            return value / counts;    
        }
        else
        {
            return 0.0;
        }
    }
    void update(const double v)
    {
        value += v;
        counts += 1.0;
    }
};

struct RidgeEnergyObject
{
    long long total_steps = 0;
    StreamingMedian streaming_median;
    StreamingMean streaming_mean;

    // Last energies that were under the threshold
    double last_energy = 0.0;
    double current_ridge = 0.0;

    // We must handle the case when the tracer STARTS above the threshold. In
    // this situation, we should not log the first time it drops below.
    bool exited_first_basin = false;

    double threshold;
    bool threshold_valid = true;

    // Keep track of the observables!
    std::vector<double> vec_means;
    std::vector<double> vec_medians;
    std::vector<double> vec_total_steps;
};


class OnePointObservables
{
protected:

    // Base objects we need
    const utils::SimulationParameters params;
    std::vector<long long> grid;
    int grid_length;
    const SpinSystem* spin_system_ptr;

    // The pointer to the last-updated point on the grid
    // Note that this is not actually a pointer =]
    unsigned int pointer = 0;

    // Define the ridge energy objects
    RidgeEnergyObject ridge_E_object;
    RidgeEnergyObject ridge_S_object;

    // Vector to track the energy
    std::vector<double> vec_energy;
    std::vector<double> vec_energy_IS;
    std::vector<std::string> vec_cache_size;
    std::vector<double> vec_acceptance_rate;
    std::vector<double> vec_walltime_per_waiting_time;

    // Other private methods. These are for ridge energies
    RidgeEnergyObject* _get_RidgeEnergyObject_ptr(const std::string which_ridge);
    void _step_ridge(const double waiting_time, const double simulation_clock, const std::string which_ridge);
    void _log_ridge(const std::string which_ridge);

public:

    // Constructor: reads in the grid from the specified grid directory
    OnePointObservables(const utils::SimulationParameters params, const SpinSystem& spin_system);

    void step(const double waiting_time, const double simulation_clock);
    json as_json() const;
};


#endif
