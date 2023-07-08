#ifndef SPIN_H
#define SPIN_H

#include <random>

#include "utils.h"
#include "energy_mapping.h"


class SpinSystem
{

// Accessible only from within the class or it children
protected:
    parameters::SimulationParameters params;
    EnergyMapping* emap_ptr;
    ap_uint<PRECISON> current_state;
    mutable parameters::SimulationStatistics sim_stats;

    // Gillespie only //////////////////////////////////////////////////////
    // Pointer to the delta E and exit rates
    double* _exit_rates = 0;
    ap_uint<PRECISON>* _neighbors = 0;
    double* _neighboring_energies = 0;
    std::vector<double> _normalized_exit_rates;
    std::exponential_distribution<double> total_exit_rate_dist;

    // Fills the exit_rates and delta_E arrays and returns the total exit
    // rate.
    double _calculate_exit_rates(const double current_energy) const;
    ////////////////////////////////////////////////////////////////////////

    // Standard only
    std::uniform_real_distribution<> uniform_0_1_distribution;
    std::uniform_int_distribution<> spin_distribution;
    ////////////////////////////////////////////////////////////////////////

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Pointer to the neighbors and neighboring energies, used in the inherent
    // structure computation
    // ap_uint<PRECISON>* neighbors = 0;
    // double* neighboring_energies = 0;


    // Initialize some objects for storing the previous and current values of
    // things. This is required for some of the 2-point observables.
    parameters::StateProperties _prev, _curr;

    /**
     * @brief [brief description]
     * @details [long description]
     */
    void _first_time_state_initialization_();

    void _init_gillespie();
    void _teardown_gillespie();
    void _init_standard();
    void _teardown_standard();

    // void _helper_fill_neighboring_energies(int *, int, double *) const;
    // long long _help_get_inherent_structure() const;
    // long long _get_inherent_structure() const;

    // // Updater for the previous values; this should be done at the end of
    // // every recording phase
    void _init_previous_state_();
    void _init_current_state_();

// Accessible outside of the class instance
public:

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param s [description]
     * @param g [description]
     */
    SpinSystem(const parameters::SimulationParameters params, EnergyMapping& emap);

    /**
     * @brief [brief description]
     * @details [long description]
     * @return [description]
     */
    double energy() const;

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param state [description]
     */
    void set_state(ap_uint<PRECISON> state);

    /**
     * @brief [brief description]
     * @details [long description]
     * @return [description]
     */
    std::string binary_state() const;

    /**
     * @brief Gets the inherent structure and logs timing information
     * @details [long description]
     * @return [description]
     */
    ap_uint<PRECISON> get_inherent_structure_();

    parameters::StateProperties get_previous_state() const {return _prev;}
    parameters::StateProperties get_current_state() const {return _curr;}
    EnergyMapping* get_emap_ptr() const {return emap_ptr;}
    parameters::SimulationStatistics get_sim_stats() const {return sim_stats;}
    // double get_average_neighboring_energy() const;
    
    double _step_standard();
    double _step_gillespie();
    double step();
    void summarize();

    ~SpinSystem();
};

#endif
