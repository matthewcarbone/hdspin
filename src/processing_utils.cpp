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


json statistics(const std::vector<json> results, const std::string key)
{
    json j;
    const int N = results.size();
    const int M = results[0].size();

    // Load in
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


namespace processing_utils
{

void postprocess()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != MASTER){return;}

    std::vector<json> results = read_all_results();

    json j;
    j["energy"] = statistics(results, "energy");
    utils::json_to_file_no_format(j, "final.json");

}

}
