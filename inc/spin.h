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
    Counter counter;

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Pointer to the neighbors and neighboring energies, used in the inherent
    // structure computation
    // ap_uint<PRECISON>* neighbors = 0;
    // double* neighboring_energies = 0;


    // Initialize some objects for storing the previous and current values of
    // things:
    parameters::StateProperties _prev, _curr;

    /**
     * @brief [brief description]
     * @details [long description]
     */
    void _first_time_state_initialization_();

    /**
     * @brief [brief description]
     * @details [long description]
     */
    void _fill_neighbors_();

    /**
     * @brief [brief description]
     * @details [long description]
     */
    void _fill_neighboring_energies_();

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

    parameters::StateProperties get_previous_state() const {return _prev;}
    parameters::StateProperties get_current_state() const {return _curr;}
    // double get_average_neighboring_energy() const;
    // ~SpinSystem();
};


class GillespieSpinSystem : public SpinSystem
{
private:

    // Pointer to the delta E and exit rates
    double *_exit_rates = 0;
    ap_uint<PRECISON> *_neighbors = 0;
    double *_neighboring_energies = 0;

    std::vector<double> _normalized_exit_rates;
    std::exponential_distribution<long double> total_exit_rate_dist;

    // Fills the exit_rates and delta_E arrays and returns the total exit
    // rate.
    double _calculate_exit_rates(const double current_energy) const;

public:
    GillespieSpinSystem(const parameters::SimulationParameters params, EnergyMapping& emap);

    // Step computes the neighboring energies, delta E values and exit rates,
    // then based on that information, steps the spin configuration and
    // returns the waiting time. Note that a Gillespie step is always accepted.
    long double step();

    ~GillespieSpinSystem();
};


class StandardSpinSystem : public SpinSystem
{
private:
    std::uniform_real_distribution<> uniform_0_1_distribution;
    std::uniform_int_distribution<> spin_distribution;

public:
    StandardSpinSystem(const parameters::SimulationParameters params, EnergyMapping& emap);

    // Step executes a possible alteration in the state, but not always. Thus,
    // the standard step actually returns whether or not the new state was
    // accepted: if there was a rejection, return false, else, if the proposed
    // state was accepted, return true.
    long double step();
    void summarize();
};




#endif
