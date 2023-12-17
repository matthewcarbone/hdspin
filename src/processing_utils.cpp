#include <fstream>
#include <filesystem>
#include <mpi.h>

#include "utils.h"
#include "processing_utils.h"

namespace fs = std::filesystem;


#define MASTER 0


json read_json(const std::string fname)
{
    std::ifstream f(fname);
    return json::parse(f);
}


std::vector<json> read_all_results()
{
    std::vector<json> all_json;
    try
    {
        for (const auto& entry : fs::directory_iterator(DATA_PATH))
        {
            if (entry.path().extension() == ".json")
            {
                all_json.push_back(read_json(entry.path().string()));
            }
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error while parsing json results: " << e.what() << std::endl;
    }
    return all_json;
}


json get_standard_statistics(const std::vector<json> results, const std::string key)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    json j;
    const size_t N = results.size();
    const size_t M = results[0][key].size();

    std::vector<std::vector<double>> target;
    for (auto &result : results){target.push_back(result[key]);}

    // For every time point, get the variables of interest
    std::vector<double> means, medians, standard_deviations, standard_errors, tmp;
    for (int jj=0; jj<M; jj++)
    {
        for (int ii=0; ii<N; ii++){tmp.push_back(target[ii][jj]);}
        const double mean = utils::mean_vector(tmp);
        const double median = utils::median_vector(tmp);
        const double variance = utils::variance_vector(tmp);
        const double standard_deviation = sqrt(variance);
        const double standard_error = standard_deviation / sqrt(N);
        means.push_back(mean);
        medians.push_back(median);
        standard_deviations.push_back(standard_deviation);
        standard_errors.push_back(standard_error);
        tmp.clear();
    }

    j["mean"] = means;
    j["median"] = medians;
    j["standard_deviation"] = standard_deviations;
    j["standard_error"] = standard_errors;

    const double dt = utils::get_time_delta(t_start);

    printf("Standard statistics : %s : done in %.02f s\n", key.c_str(), dt);

    return j;
}


json get_ridge_statistics(const std::vector<json> results, const std::string key)
{

    auto t_start = std::chrono::high_resolution_clock::now();

    const std::string mean_key = key + "_mean";
    const std::string median_key = key + "_median";
    const std::string total_steps_key = key + "_total_steps";

    json j;
    const size_t N = results.size();
    const size_t M = results[0][mean_key].size();

    // std::vector<std::vector<double>> target = _parse_results(results, key);
    // std::vector<double> means, medians, standard_deviations, standard_errors;
    // std::vector<double> tmp_means, tmp_medians, tmp_weights;

    // Load in the ridge energy statistics, we'll focus on E first
    std::vector<std::vector<double>> ridge_means, ridge_medians, ridge_weights;
    
    // std::cout << mean_key << " " << median_key << " " << total_steps_key << std::endl;
    for (auto &result : results)
    {
        ridge_means.push_back(result[mean_key]);
        ridge_medians.push_back(result[median_key]);
        ridge_weights.push_back(result[total_steps_key]);
    }

    // Logic here is a bit different than standard statistics, as we need to 
    // compute weighted averages and weighted variances
    // this information is saved:
    //
    // j["ridge_E_mean"] = ridge_E_object.vec_means;
    // j["ridge_E_median"] = ridge_E_object.vec_medians;
    // j["ridge_E_total_steps"] = ridge_E_object.vec_total_steps;

    std::vector<double> mean_avg, mean_standard_deviations, mean_standard_errors;
    std::vector<double> median_avg, median_standard_deviations, median_standard_errors;
    std::vector<double> all_weights;

    std::vector<double> tmp_ridge_means, tmp_ridge_medians, tmp_ridge_weights;
    std::vector<double> tmp_var;
    double total_N;
    for (int jj=0; jj<M; jj++)
    {

        // Fill all the tracers
        for (int ii=0; ii<N; ii++)
        {   

            // Special logic here. We do not want to count anything in which
            // the value for the ridge is 0 and the weight is zero. This is a
            // case in which nothing has been logged yet. The equivalent in the
            // old python postprocessing script was to use a masked array.
            // Here we have to handle it explicitly

            const double _mean = ridge_means[ii][jj];
            const double _median = ridge_medians[ii][jj];
            const double _weight = ridge_weights[ii][jj];

            if ((_mean == 0.0) && (_median == 0.0) && (_weight == 0.0))
            {
                continue;
            }

            tmp_ridge_means.push_back(_mean);
            tmp_ridge_medians.push_back(_median);
            tmp_ridge_weights.push_back(_weight);
        }

        // std::cout << "tmp_ridge_means.size()==" << tmp_ridge_means.size() << std::endl;
        // std::cout << "tmp_ridge_medians.size()==" << tmp_ridge_medians.size() << std::endl;
        // std::cout << "tmp_ridge_weights.size()==" << tmp_ridge_weights.size() << std::endl;

        // Special case here
        if (tmp_ridge_weights.size() == 0)
        {
            all_weights.push_back(0.0);
            mean_avg.push_back(0.0);
            mean_standard_deviations.push_back(0.0);
            mean_standard_errors.push_back(0.0);
            median_avg.push_back(0.0);
            median_standard_deviations.push_back(0.0);
            median_standard_errors.push_back(0.0);
        }

        else
        {
            total_N = std::accumulate(tmp_ridge_weights.begin(), tmp_ridge_weights.end(), 0.0);
            all_weights.push_back(total_N);
            // Do logic for the mean ridge energies here
            double mu = utils::weighted_mean_vector(tmp_ridge_means, tmp_ridge_weights);

            // Calculate the variance of the means
            for (int kk=0; kk<tmp_ridge_means.size(); kk++)  // problem is here I think
                // Need to account for the fact that tmp_ridge_means can have
                // sizes different from N
            {
                tmp_var.push_back(pow(mu - tmp_ridge_means[kk], 2));
            }
            double var = utils::weighted_mean_vector(tmp_var, tmp_ridge_weights);
            double std = sqrt(var);
            
            mean_avg.push_back(mu);
            mean_standard_deviations.push_back(std);
            mean_standard_errors.push_back(std / sqrt(total_N));
            tmp_var.clear();


            // Do logic for the median ridge energies here
            mu = utils::weighted_mean_vector(tmp_ridge_medians, tmp_ridge_weights);

            // Calculate the variance of the medians
            for (int kk=0; kk<tmp_ridge_medians.size(); kk++)
            {
                tmp_var.push_back(pow(mu - tmp_ridge_medians[kk], 2));
            }
            var = utils::weighted_mean_vector(tmp_var, tmp_ridge_weights);
            std = sqrt(var);
            double stderr = std / sqrt(total_N);
            
            median_avg.push_back(mu);
            median_standard_deviations.push_back(std);
            median_standard_errors.push_back(stderr);
            tmp_var.clear();
        }

        tmp_ridge_means.clear();
        tmp_ridge_medians.clear();
        tmp_ridge_weights.clear();
    }

    j["mean_mean"] = mean_avg;
    j["mean_std"] = mean_standard_deviations;
    j["mean_std_err"] = mean_standard_errors;
    j["all_weights"] = all_weights;
    j["median_mean"] = median_avg;
    j["median_std"] = median_standard_deviations;
    j["median_std_err"] = median_standard_errors;

    const double dt = utils::get_time_delta(t_start);

    printf("Ridge statistics : %s : done in %.02f s\n", key.c_str(), dt);

    return j;
}


json get_aging_config_statistics(const std::vector<json> results)
{

    auto t_start = std::chrono::high_resolution_clock::now();

    json j;
    const size_t N = results.size();
    const double sqrt_N = sqrt(N);
    const size_t M = results[0]["aging_config_pi1"].size();

    std::vector<std::vector<std::string>> aging_pi1, aging_pi2;
    for (auto &result : results)
    {
        aging_pi1.push_back(result["aging_config_pi1"]);
        aging_pi2.push_back(result["aging_config_pi2"]);
    }

    std::vector<double> tmp;
    std::vector<double> means, medians, standard_deviations, standard_errors;
    for (size_t jj=0; jj<M; jj++)
    {
        tmp.clear();
        for (size_t ii=0; ii<N; ii++)
        {
            tmp.push_back(aging_pi1[ii][jj] == aging_pi2[ii][jj] ? 1.0 : 0.0);
        }
        const double mean = utils::mean_vector(tmp);
        const double median = utils::median_vector(tmp);
        const double variance = utils::variance_vector(tmp);
        means.push_back(mean);
        medians.push_back(median);
        const double standard_deviation = sqrt(variance);
        const double standard_error = standard_deviation / sqrt_N;
        standard_deviations.push_back(standard_deviation);
        standard_errors.push_back(standard_error);
    }

    j["mean"] = means;
    j["median"] = medians;
    j["standard_deviation"] = standard_deviations;
    j["standard_error"] = standard_errors;

    const double dt = utils::get_time_delta(t_start);

    printf("Aging config statistics : done in %.02f s\n", dt);

    return j;
}



json get_aging_basin_statistics(const std::vector<json> results, const std::string key)
{
    // key == aging_basin_E or aging_basin_S
    auto t_start = std::chrono::high_resolution_clock::now();

    const std::string key1 = key + "_index_1";
    const std::string key2 = key + "_index_2";
    const std::string key_in = key + "_prev_state_in_basin_1";

    json j;
    const size_t N = results.size();
    const size_t M = results[0][key1].size();

    std::vector<std::vector<int>> aging_pi1, aging_pi2;
    std::vector<std::vector<int>> in_basin;
    for (auto &result : results)
    {
        aging_pi1.push_back(result[key1]);
        aging_pi2.push_back(result[key2]);
        in_basin.push_back(result[key_in]);
    }

    std::vector<double> tmp;
    std::vector<double> means, medians, standard_deviations, standard_errors;
    for (size_t jj=0; jj<M; jj++)
    {
        tmp.clear();
        for (size_t ii=0; ii<N; ii++)
        {
            if (in_basin[ii][jj] == 1)
            {
                tmp.push_back(aging_pi1[ii][jj] == aging_pi2[ii][jj] ? 1.0 : 0.0);
            }
            
        }

        if (tmp.size() == 0)
        {
            means.push_back(0.0);
            medians.push_back(0.0);
            standard_deviations.push_back(0.0);
            standard_errors.push_back(0.0);
        }
        else
        {
            const double mean = utils::mean_vector(tmp);
            const double median = utils::median_vector(tmp);
            const double variance = utils::variance_vector(tmp);
            means.push_back(mean);
            medians.push_back(median);
            const double standard_deviation = sqrt(variance);
            const double standard_error = standard_deviation / sqrt(tmp.size());
            standard_deviations.push_back(standard_deviation);
            standard_errors.push_back(standard_error);
        }
    }

    j["mean"] = means;
    j["median"] = medians;
    j["standard_deviation"] = standard_deviations;
    j["standard_error"] = standard_errors;

    const double dt = utils::get_time_delta(t_start);

    printf("Aging basin statistics : done in %.02f s\n", dt);

    return j;
}    



json load_grids()
{
    json j;

    std::vector<double> grid_energy;
    utils::load_grid_(grid_energy, ENERGY_GRID_PATH);
    j["energy"] = grid_energy;

    std::vector<double> grid_pi1;
    utils::load_grid_(grid_pi1, PI1_GRID_PATH);
    j["pi1"] = grid_pi1;

    std::vector<double> grid_pi2;
    utils::load_grid_(grid_pi2, PI2_GRID_PATH);
    j["pi2"] = grid_pi2;

    return j;
}


std::vector<double> get_psi_grid(const size_t N)
{
    const double e = 2.71828182845904523536;
    std::vector<double> v;
    for (size_t ii=0; ii<N; ii++){v.push_back(pow(e, ii));}
    return v;
}


namespace processing_utils
{

void postprocess()
{

    MPI_Barrier(MPI_COMM_WORLD);  // <----- Barrier to make sure all ranks finish

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != MASTER){return;}

    const json config = read_json(CONFIG_PATH);
    const bool valid_entropic_attractor = config["valid_entropic_attractor"];

    printf("Postprocessing results\n");

    const std::vector<json> results = read_all_results();

    json j;
    j["energy"] = get_standard_statistics(results, "energy");
    j["acceptance_rate"] = get_standard_statistics(results, "acceptance_rate");
    j["cache_size"] = get_standard_statistics(results, "cache_size");
    j["walltime_per_waitingtime"] = get_standard_statistics(results, "walltime_per_waitingtime");
    j["emax"] = get_standard_statistics(results, "emax");
    j["ridge_E"] = get_ridge_statistics(results, "ridge_E");
    if (valid_entropic_attractor)
    {
        j["ridge_S"] = get_ridge_statistics(results, "ridge_S");
    }

    j["psi_config"] = get_standard_statistics(results, "psi_config");
    j["psi_basin_E"] = get_standard_statistics(results, "psi_basin_E");
    if (valid_entropic_attractor)
    {
        j["psi_basin_S"] = get_standard_statistics(results, "psi_basin_S");
    }

    j["aging_config"] = get_aging_config_statistics(results);
    j["aging_basin_E"] = get_aging_basin_statistics(results, "aging_basin_E");
    if (valid_entropic_attractor)
    {
        j["aging_basin_S"] = get_aging_basin_statistics(results, "aging_basin_S");
    }

    j["grids"] = load_grids();
    j["grids"]["psi"] = get_psi_grid(j["psi_config"]["mean"].size());
    

    utils::json_to_file(j, RESULTS_PATH);
    printf("Results saved to %s\n", RESULTS_PATH);
    utils::cleanup_directories();
    printf("Simulation completed successfully\n");
}

}
