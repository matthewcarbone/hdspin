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
    json j;
    const int N = results.size();
    const int M = results[0][key].size();

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
    return j;
}


json get_ridge_statistics(const std::vector<json> results, const std::string key)
{

    const std::string mean_key = key + "_mean";
    const std::string median_key = key + "_median";
    const std::string total_steps_key = key + "_total_steps";

    json j;
    const int N = results.size();
    const int M = results[0][mean_key].size();

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
            tmp_ridge_means.push_back(ridge_means[ii][jj]);
            tmp_ridge_medians.push_back(ridge_medians[ii][jj]);
            tmp_ridge_weights.push_back(ridge_weights[ii][jj]);
        }

        total_N = std::accumulate(tmp_ridge_weights.begin(), tmp_ridge_weights.end(), 0.0);
        all_weights.push_back(total_N);


        // Do logic for the mean ridge energies here

        double mu = utils::weighted_mean_vector(tmp_ridge_means, tmp_ridge_weights);

        // Calculate the variance of the means
        
        for (int ii=0; ii<N; ii++)
        {
            tmp_var.push_back(pow(mu - ridge_means[ii][jj], 2));
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
        for (int ii=0; ii<N; ii++)
        {
            tmp_var.push_back(pow(mu - ridge_medians[ii][jj], 2));
        }
        var = utils::weighted_mean_vector(tmp_var, tmp_ridge_weights);
        std = sqrt(var);
        double stderr = std / sqrt(total_N);
        
        median_avg.push_back(mu);
        median_standard_deviations.push_back(std);
        median_standard_errors.push_back(stderr);
        tmp_var.clear();

        tmp_ridge_means.clear();
        tmp_ridge_medians.clear();
        tmp_ridge_weights.clear();
    }

    std::cout << "made it to the end " << std::endl;

    j["mean_mean"] = mean_avg;
    j["mean_std"] = mean_standard_deviations;
    j["mean_std_err"] = mean_standard_errors;
    j["all_weights"] = all_weights;
    j["median_mean"] = median_avg;
    j["median_std"] = median_standard_deviations;
    j["median_std_err"] = median_standard_errors;

    return j;
}


namespace processing_utils
{

void postprocess()
{

    MPI_Barrier(MPI_COMM_WORLD);  // <----- Barrier to make sure all finish

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != MASTER){return;}

    std::vector<json> results = read_all_results();

    json j;
    j["energy"] = get_standard_statistics(results, "energy");
    j["acceptance_rate"] = get_standard_statistics(results, "acceptance_rate");
    j["walltime_per_waitingtime"] = get_standard_statistics(results, "walltime_per_waitingtime");
    j["emax"] = get_standard_statistics(results, "emax");
    j["ridge_E"] = get_ridge_statistics(results, "ridge_E");

    utils::json_to_file(j, "final.json");
    std::cout << "Averaged results saved to final.json" << std::endl;
}

}
