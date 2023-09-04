#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <vector>
#include <string>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix read_file(const std::string& filename);
std::vector<std::string> get_all_results_filenames(const std::string& directory);
void save_matrices_to_file(const std::string& save_path, const Matrix& matrix);
double calculate_mean(const Vector& vec);
double calculate_std_dev(const Vector& vec, double mean);
double calculate_sem(double std_dev, size_t n);
double calculate_weighted_mean(const Vector& values, const Vector& weights);
double calculate_weighted_variance(const Vector& values, const Vector& weights, double weighted_mean);
double calculate_weighted_sem(double weighted_var, const Vector& weights);
void obs1(const std::vector<std::string>& all_filenames, const std::string& substring, const std::string& save_path);
void ridge(const std::vector<std::string>& all_filenames, const std::string& substring, const std::string& save_path);
void cache_size(const std::vector<std::string>& all_filenames, const std::string& substring, const std::string& save_path);

#endif
