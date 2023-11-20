#ifndef TEST_OBS1_H
#define TEST_OBS1_H

#include <algorithm>
#include <random>
// #include <stdio.h>

#include "obs1.h"


namespace test_obs1
{

bool test_streaming_median()
{

    std::default_random_engine generator;
    generator.seed(123);
    std::uniform_real_distribution<double> distribution(-10.0, 10.0);

    std::vector<double> v;
    double true_median;

    for (unsigned int vector_size=2; vector_size<10000; vector_size+=101)
    {

        v.clear();
        
        // Fill the vector
        for (unsigned int ii=0; ii<vector_size; ii++)
        {
            const double val = distribution(generator);
            v.push_back(val);
        }

        // Compute the median via sorting first
        std::sort(v.begin(), v.end());

        const unsigned int tmp_ii = vector_size / 2;

        // Even
        // 0 1 2 3 4 5 need index 2 and 3 for size 6
        // 0 1 2 3 4 5 6 7 need index 3 and 4 for size 8
        // 0 1 2 ... N need index N/2 and N/2-1 for size N
        // Where N is even
        if (vector_size % 2 == 0)
        {
            true_median = (v[tmp_ii] + v[tmp_ii-1]) / 2.0;
        }

        // Odd
        // 0 1 2 need index 1 for size 3
        // 0 1 2 3 4 need index 2 for size 5
        // 0 1 2 ... N need index N/2 (integer division) for size N
        // Where N is odd
        else
        {
            true_median = v[tmp_ii];
        }

        // Now get the median calculated using the rolling method
        std::shuffle(std::begin(v), std::end(v), generator);
        StreamingMedian streaming_median = StreamingMedian();
        for (unsigned int ii=0; ii<vector_size; ii++)
        {
            streaming_median.update(v[ii]); // Update in random order
        }    

        const double calc_median = streaming_median.median();

        // printf("%i: True vs. Calc %.04f vs. %.04f\n", vector_size, true_median, calc_median);
        if (calc_median != true_median){return false;}
    }

    return true;
}

}

#endif
