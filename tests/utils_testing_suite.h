#ifndef UTILS_TESTING_SUITE_H
#define UTILS_TESTING_SUITE_H

#include <vector>
#include <numeric>

double _mean_vector(const std::vector<double> v)
{
    const double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}

double _variance_vector(const std::vector<double> v)
{
    const double mean = _mean_vector(v);
    const double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    return sq_sum / v.size() - mean * mean;
}

#endif
