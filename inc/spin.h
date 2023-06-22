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

    // Initialize the MT random number generator and seed with random_device
    // This is seeded in the constructor
    mutable std::mt19937 generator;

    // Pointer to the neighbors and neighboring energies, used in the inherent
    // structure computation
    ap_uint<PRECISON>* neighbors = 0;
    double* neighboring_energies = 0;


    // Initialize some objects for storing the previous and current values of
    // things:
    // Vals prev, curr;

    // Multiplier for sampling from the waiting time and total
    // exit rate. This is 1 by default but is set to try and find the
    // equivalent Gillespie simulation for the "loop" standard dynamics.
    double _waiting_time_multiplier = 1.0;

    /**
     * @brief [brief description]
     * @details [long description]
     */
    void _initialize_state_();

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
    // void _init_prev();
    // void _init_curr();

// Accessible outside of the class instance
public:

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param s [description]
     * @param g [description]
     */
    SpinSystem(const parameters::SimulationParameters, EnergyMapping&);

    /**
     * @brief [brief description]
     * @details [long description]
     * @return [description]
     */
    double energy() const;

    /**
     * @brief [brief description]
     * @details [long description]
     * @return [description]
     */
    std::string binary_state() const;

    // Vals get_prev() const {return prev;}
    // Vals get_curr() const {return curr;}
    // double get_average_neighboring_energy() const;
    ~SpinSystem();
};


// class GillespieSpinSystem : public SpinSystem
// {
// private:

//     // Pointer to the delta E and exit rates
//     double *_delta_E = 0;
//     double *_exit_rates = 0;

//     std::vector<double> _normalized_exit_rates;
//     std::exponential_distribution<long double> total_exit_rate_dist;

//     // Fills the exit_rates and delta_E arrays and returns the total exit
//     // rate.
//     double _calculate_exit_rates() const;

// public:
//     GillespieSpinSystem(const RuntimeParameters, EnergyMapping&);

//     // Step computes the neighboring energies, delta E values and exit rates,
//     // then based on that information, steps the spin configuration and
//     // returns the waiting time. Note that a Gillespie step is always accepted.
//     long double step();

//     ~GillespieSpinSystem();
// };


// class StandardSpinSystem : public SpinSystem
// {
// private:
//     std::uniform_real_distribution<> uniform_0_1_distribution;
//     std::uniform_int_distribution<> spin_distribution;

// public:
//     StandardSpinSystem(const RuntimeParameters, EnergyMapping&);

//     // Step executes a possible alteration in the state, but not always. Thus,
//     // the standard step actually returns whether or not the new state was
//     // accepted: if there was a rejection, return false, else, if the proposed
//     // state was accepted, return true.
//     long double step();
// };




#endif
