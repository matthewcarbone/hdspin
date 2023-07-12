#ifndef OBS1_H
#define OBS1_H

#include <queue>
#include <vector>
#include <unordered_set>

#include "spin.h"
#include "utils.h"


class RollingMedian
{
protected:

    // Standard  max priority queue
    std::priority_queue <double> max_heap;

    // Reverse (min) priority queue
    std::priority_queue <double, std::vector <double>, std::greater<double>> min_heap;

public:
    RollingMedian();
    double median() const;
    void update(const double v);
};


class ObsBase
{
protected:
    const parameters::FileNames fnames;
    const parameters::SimulationParameters params;
    std::vector<long long> grid;
    int grid_length;
    const SpinSystem* spin_system_ptr;

    // The pointer to the last-updated point on the grid
    unsigned int pointer = 0;

public:
    ObsBase(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);
};

class RidgeBase : public ObsBase
{
protected:

    FILE* outfile;

    long double _ridge_energy_accumulator = 0.0;
    long long _total_steps = 0;

    // Last energies that were under the threshold
    double _last_energy = 0.0;
    double _current_ridge = 0.0;

    // We must handle the case when the tracer STARTS above the threshold. In
    // this situation, we should not log the first time it drops below.
    bool _exited_first_basin = false;

    double _threshold;
    bool _threshold_valid = true;

public:

    // Constructor: reads in the grid from the specified grid directory
    RidgeBase(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);

    // Step the grid by performing the following steps:
    // 1) Stepping the pointer
    // 2) Saving the configuration/energy information to disk
    void step(const double waiting_time, const double simulation_clock);

    ~RidgeBase();
};



class RidgeE : public RidgeBase
{
public:
    RidgeE(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);
};

class RidgeS : public RidgeBase
{
public:
    RidgeS(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);
};



class OnePointObservables : public ObsBase
{
protected:

    // Output files and pointers
    FILE* outfile_energy;
    FILE* outfile_energy_IS;
    FILE* outfile_capacity;
    FILE* outfile_acceptance_rate;
    FILE* outfile_inherent_structure_timings;
    FILE* outfile_walltime_per_waitingtime;

public:

    // Constructor: reads in the grid from the specified grid directory
    OnePointObservables(const parameters::FileNames fnames, const parameters::SimulationParameters params, const SpinSystem& spin_system);

    void step(const double waiting_time, const double simulation_clock);
    ~OnePointObservables();
};


#endif
