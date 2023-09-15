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

struct _RidgeEnergyObjects
{
    FILE* outfile;

    // long double _ridge_energy_accumulator = 0.0;
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
};


class OnePointObservables
{
protected:

    // Base objects we need
    const parameters::FileNames fnames;
    const parameters::SimulationParameters params;
    std::vector<long long> grid;
    int grid_length;
    const SpinSystem* spin_system_ptr;

    // The pointer to the last-updated point on the grid
    unsigned int pointer = 0;

    // Output files and pointers
    FILE* outfile_energy;
    FILE* outfile_energy_IS;
    FILE* outfile_capacity;
    FILE* outfile_acceptance_rate;
    FILE* outfile_inherent_structure_timings;
    FILE* outfile_walltime_per_waitingtime;

    // Define the ridge energy objects
    _RidgeEnergyObjects ridge_E_objects;
    _RidgeEnergyObjects ridge_S_objects;

    // Other private methods. These are for ridge energies
    _RidgeEnergyObjects* _get_ridge_pointer(const std::string which_ridge);
    void _step_ridge(const double waiting_time, const double simulation_clock, const std::string which_ridge);
    void _ridge_writeout(const std::string which_ridge);

public:

    // Constructor: reads in the grid from the specified grid directory
    OnePointObservables(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);

    void step(const double waiting_time, const double simulation_clock);
    ~OnePointObservables();
};


#endif
