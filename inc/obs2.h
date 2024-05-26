#ifndef OBS_PSI_H
#define OBS_PSI_H

#include <vector>
#include <set>

#include "spin.h"
#include "utils.h"


class PsiConfig
{
protected:

    const utils::SimulationParameters params;
    const SpinSystem* spin_system_ptr;
    
    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    size_t _max_counter;
    std::vector<long long> _counter;

    // Keep track internally of the waiting time for both the standard
    // trajectory and the inherent structure
    double _waiting_time = 0.0;

    size_t _out_of_counter = 0;

public:

    PsiConfig(const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double current_waiting_time);
    json as_json() const;

};

struct PsiBasinData
{
    std::vector<long long> _counter;
    std::vector<long long> _counter_unique_configs_per_basin;
    std::set<std::string> _tmp_unique_configs_in_basin;
    
    double threshold;
    bool threshold_valid = true;

    double _waiting_time = 0.0;
};

class PsiBasin
{
protected:

    // TODO- refeactor so these inherit from a base class
    const utils::FileNames fnames;
    const utils::SimulationParameters params;
    const SpinSystem* spin_system_ptr;
    size_t _max_counter;
    PsiBasinData data_E, data_S;

    void _init_E_data();
    void _init_S_data();
    // void _init_data_objects(PsiBasinData &data);
    PsiBasinData* _get_PsiBasinData_ptr(const std::string which);
    void _step(const double current_waiting_time, const std::string which);
    

public:
    PsiBasin(const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double current_waiting_time);
    json as_json() const;
};


class Aging
{
protected:
    const SpinSystem* spin_system_ptr;

    const utils::SimulationParameters params;

    // The pi 1 and 2 grids
    std::vector<double> grid_pi1;
    std::vector<double> grid_pi2;
    size_t length;

public:
    Aging(const utils::SimulationParameters params, const SpinSystem& spin_system);
};


class AgingConfigData
{
public:
    // Results holder
    std::vector<std::string> results;

    // Define the pointer
    size_t pointer = 0;


};

class AgingConfig : public Aging
{
protected:

    AgingConfigData d1, d2;
    unsigned long long state_index = 0;

    // Helpers
    void _help_step(const double simulation_clock, const size_t index);

public:
    AgingConfig(const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double simulation_clock);
    json as_json() const;
};


struct AgingBasinData
{
    // Define the pointers
    size_t pointer1 = 0;
    size_t pointer2 = 0;

    std::vector<long long> vec_basin_index_1;
    std::vector<long long> vec_basin_index_2;
    std::vector<size_t> vec_prev_state_in_basin_1;
    std::vector<size_t> vec_prev_state_in_basin_2;

    double threshold;
    bool threshold_valid = true;

    // The basin indexes reference the last basin that a tracer was in, or the
    // basin the tracer is currently in. Whether or not a tracer is currently
    // in a basin or not can be determined by comparing the energies to the
    // thresholds in the RuntimeParameters, rtp.
    long long basin_index_1 = 0;  // First basin index
    long long basin_index_2 = 0;  // Second basin index
};


class AgingBasin : public Aging
{
protected:
    AgingBasinData data_E, data_S;

    void _init_E_data();
    void _init_S_data();
    AgingBasinData* _get_PsiBasinData_ptr(const std::string which);
    void _help_step_1_(const double simulation_clock, const double prev_energy, AgingBasinData* data_ptr);
    void _help_step_2_(const double simulation_clock, const double prev_energy, AgingBasinData* data_ptr);
    void _help_step(const double simulation_clock, const std::string which);

public:
    AgingBasin(const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double simulation_clock);
    json as_json();
};

#endif
