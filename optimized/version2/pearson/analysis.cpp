/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <list>
#include <vector>

namespace Analysis {

std::vector<double> correlation_coefficients(std::vector<Vector> &datasets) {
    std::vector<double> result{};
    int datasetSize = datasets.size();
    std::vector<double> means(datasetSize);
    std::vector<double> mm_over(datasetSize);
    // Pre process the mean values. Previously it calculated for each
    // combination.
    for (auto i{0}; i < datasetSize; i++) {
        auto data = datasets[i].get_data();
        auto size = datasets[i].get_size();
        means[i] = datasets[i].mean();
        double sum = 0.0;
        for (auto j{0}; j < size; j++) {
            double deriviation = data[j] - means[i];
            sum += deriviation * deriviation;
        }
        mm_over[i] = sum;
    }

    for (auto sample1{0}; sample1 < datasetSize - 1; sample1++) {
        for (auto sample2{sample1 + 1}; sample2 < datasetSize; sample2++) {
            auto corr{pearson(datasets[sample1], datasets[sample2],
                              means[sample1], means[sample2], mm_over[sample1],
                              mm_over[sample2])};
            result.push_back(corr);
        }
    }
    return result;
}

double pearson(Vector &vec1, Vector &vec2, double x_mean, double y_mean,
               double x_mm_over_x, double y_mm_over_y) {

    // Pre fetch data to remove overide indexing calls.
    auto x = vec1.get_data();
    auto y = vec2.get_data();
    double convSum = 0.0;
    int size = vec1.get_size();

    for (auto i{0}; i < size; i++) {
        convSum += (x[i] - x_mean) * (y[i] - y_mean);
    }

    const double denom = std::sqrt(x_mm_over_x * y_mm_over_y);

    const double r = convSum / denom;

    return std::max(std::min(r, 1.0), -1.0);
}

}; // namespace Analysis
