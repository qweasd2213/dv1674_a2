/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
std::vector<double> correlation_coefficients(std::vector<Vector> &datasets);
double pearson(Vector &vec1, Vector &vec2, double x_mean, double y_mean,
               double xmm_over, double ymm_over);
}; // namespace Analysis

#endif
