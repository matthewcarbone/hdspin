#ifndef OBS_PSI_H
#define OBS_PSI_H

#include <vector>
#include <set>

#include "spin.h"
#include "utils.h"


class PsiConfig
{
protected:

    const utils::FileNames fnames;
    const utils::SimulationParameters params;
    FILE *outfile;
    const SpinSystem* spin_system_ptr;
    
    // The (roughly) maximum timestep on the counter. It's padded at the end,
    // and is taken care of during post processing.
    long long _max_counter;
    std::vector<long long> _counter;

    // Keep track internally of the waiting time for both the standard
    // trajectory and the inherent structure
    double _waiting_time = 0.0;

    int _out_of_counter = 0;

public:

    PsiConfig(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double current_waiting_time);
    ~PsiConfig();

};

struct PsiBasinData
{
    FILE* outfile;

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
    long long _max_counter;
    PsiBasinData E_data, S_data;

    void _init_E_data();
    void _init_S_data();
    // void _init_data_objects(PsiBasinData &data);
    PsiBasinData* _get_psi_basin_data_pointer(const std::string which);
    void _step(const double current_waiting_time, const std::string which);
    void _dump_outfile(const std::string which);
    

public:
    PsiBasin(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double current_waiting_time);
    ~PsiBasin();
};


class Aging
{
protected:
    const SpinSystem* spin_system_ptr;

    const utils::FileNames fnames;
    const utils::SimulationParameters params;

    // The pi 1 and 2 grids
    std::vector<long long> grid_pi1;
    std::vector<long long> grid_pi2;
    int length;

public:
    Aging(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system);
};


class AgingConfig : public Aging
{
protected:

    // Results
    std::vector<std::string> results1;
    std::vector<std::string> results2;

    // Define the pointers
    int pointer1 = 0;
    int pointer2 = 0;

    // The outstream for this tracker
    FILE* outfile;

    // Helpers
    void _help_step_1(const double simulation_clock);
    void _help_step_2(const double simulation_clock);
    // void _help_step_2(const double, const long long);

    // void _dump_outfile();

public:
    AgingConfig(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double simulation_clock);
    ~AgingConfig();
};



struct AgingBasinData
{
    FILE* outfile;

    // Define the pointers
    int pointer1 = 0;
    int pointer2 = 0;

    std::vector<long long> vec_basin_index_1;
    std::vector<long long> vec_basin_index_2;
    std::vector<int> vec_prev_state_in_basin_1;
    std::vector<int> vec_prev_state_in_basin_2;

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
    AgingBasinData E_data, S_data;

    void _init_E_data();
    void _init_S_data();
    AgingBasinData* _get_psi_basin_data_pointer(const std::string which);
    void _dump_outfile(const std::string which);
    void _help_step_1_(const double simulation_clock, const double prev_energy, AgingBasinData* data_ptr);
    void _help_step_2_(const double simulation_clock, const double prev_energy, AgingBasinData* data_ptr);
    void _help_step(const double simulation_clock, const std::string which);

public:
    AgingBasin(const utils::FileNames fnames, const utils::SimulationParameters params, const SpinSystem& spin_system);
    void step(const double simulation_clock);
    ~AgingBasin();
};



#endif
