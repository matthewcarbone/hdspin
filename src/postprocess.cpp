#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <postprocess.h>

namespace fs = std::filesystem;

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Read a file into a Matrix
Matrix read_file(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Could not open the file: " << filename << std::endl;
        return Matrix();
    }
    Matrix matrix;
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        Vector row;
        double val;
        while (iss >> val) {
            row.push_back(val);
        }
        matrix.push_back(row);
    }
    return matrix;
}

// Get all filenames in a directory
std::vector<std::string> get_all_results_filenames(const std::string& directory) {
    std::vector<std::string> filenames;
    if (!fs::exists(directory) || !fs::is_directory(directory)) {
        std::cerr << "Directory does not exist or is not readable: " << directory << std::endl;
        return std::vector<std::string>();
    }
    for (const auto& entry : fs::directory_iterator(directory)) {
        filenames.push_back(entry.path().string());
    }
    return filenames;
}

// Save matrices to a file
void save_matrices_to_file(const std::string& save_path, const Matrix& matrix) {
    std::ofstream outfile(save_path);
    if (!outfile.is_open()) {
        std::cerr << "Could not open the file for writing: " << save_path << std::endl;
        return;
    }

    // Write the matrix to file in scientific notation with 8 decimal places
    outfile << std::scientific << std::setprecision(8);
    for (const auto& row : matrix) {
        for (const auto& val : row) {
            outfile << val << " ";
        }
        outfile << "\n";
    }

    // Check for write errors
    if (!outfile.good()) {
        std::cerr << "Error occurred while writing to the file: " << save_path << std::endl;
        return;
    }

    outfile.close();

    // Check for close errors
    if (outfile.fail()) {
        std::cerr << "Error occurred while closing the file: " << save_path << std::endl;
    }
}

// Calculate mean of a vector
double calculate_mean(const Vector& vec) {
    if (vec.empty()) return 0;
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

// Calculate standard deviation of a vector
double calculate_std_dev(const Vector& vec, double mean) {
    if (vec.empty()) return 0;
    double sum_sq_diff = 0.0;
    for (const auto& val : vec) {
        sum_sq_diff += std::pow(val - mean, 2);
    }
    return std::sqrt(sum_sq_diff / vec.size());
}

// Calculate standard error of the mean (SEM)
double calculate_sem(double std_dev, size_t n) {
    if (n == 0) return 0;
    return std_dev / std::sqrt(static_cast<double>(n));
}

double calculate_weighted_mean(const Vector& values, const Vector& weights) {
    double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sum_weights == 0) return std::nan("");

    double weighted_sum = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        weighted_sum += values[i] * weights[i];
    }
    return weighted_sum / sum_weights;
}

// Calculate weighted variance
double calculate_weighted_variance(const Vector& values, const Vector& weights, double weighted_mean) {
    double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sum_weights == 0) return std::nan(""); 

    double weighted_sum_sq_diff = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        weighted_sum_sq_diff += weights[i] * std::pow(values[i] - weighted_mean, 2);
    }
    return weighted_sum_sq_diff / sum_weights;
}

// Calculate weighted standard error of the mean (SEM)
double calculate_weighted_sem(double weighted_var, const Vector& weights) {
    double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sum_weights == 0) return std::nan("");

    return std::sqrt(weighted_var) / std::sqrt(sum_weights);
}

void obs1(const std::vector<std::string>& all_filenames, const std::string& substring, const std::string& save_path) {
    std::vector<Matrix> matrices;
    for (const auto& filename : all_filenames) {
        if (filename.find(substring) != std::string::npos) {
            Matrix mat = read_file(filename);
            if (mat.empty() || mat[0].empty()) {
                std::cerr << "Warning: Empty or invalid matrix read from file: " << filename << std::endl;
                continue;
            }
            matrices.push_back(mat);
        }
    }
    
    if (matrices.empty()) {
        std::cerr << "Error: No valid matrices found for substring: " << substring << std::endl;
        return;
    }
    
    size_t rows = matrices[0].size();
    size_t cols = matrices[0][0].size();

    // Check for consistent matrix dimensions
    for (const auto& matrix : matrices) {
        if (matrix.size() != rows || matrix[0].size() != cols) {
            std::cerr << "Error: Inconsistent matrix dimensions." << std::endl;
            return;
        }
    }

    // Further computation and saving logic
    Matrix mu(rows, Vector(cols, 0.0));
    Matrix sd(rows, Vector(cols, 0.0));
    Matrix stderr(rows, Vector(cols, 0.0));

    // Calculate mean, standard deviation, and standard error
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            Vector column_values;
            for (const auto& matrix : matrices) {
                column_values.push_back(matrix[i][j]);
            }
            
            // Compute mean using utility function
            double mean = calculate_mean(column_values);
            mu[i][j] = mean;

            // Compute standard deviation using utility function
            double std_dev = calculate_std_dev(column_values, mean);
            sd[i][j] = std_dev;

            // Compute standard error using utility function
            double sem = calculate_sem(std_dev, matrices.size());
            stderr[i][j] = sem;
        }
    }

    save_matrices_to_file(save_path, mu);
}

void ridge(const std::vector<std::string>& all_filenames, const std::string& substring, const std::string& save_path) {
    std::vector<Matrix> matrices;
    for (const auto& filename : all_filenames) {
        if (filename.find(substring) != std::string::npos) {
            Matrix mat = read_file(filename);
            if (mat.empty() || mat[0].empty()) {
                std::cerr << "Warning: Empty or invalid matrix read from file: " << filename << std::endl;
                continue;
            }
            matrices.push_back(mat);
        }
    }
    
    if (matrices.empty()) {
        std::cerr << "Error: No valid matrices found for substring: " << substring << std::endl;
        return;
    }
    
    size_t rows = matrices[0].size();
    size_t cols = matrices[0][0].size();

    // Check for consistent matrix dimensions
    for (const auto& matrix : matrices) {
        if (matrix.size() != rows || matrix[0].size() != cols) {
            std::cerr << "Error: Inconsistent matrix dimensions." << std::endl;
            return;
        }
    }

    Vector mu1(rows, 0.0), mu2(rows, 0.0);
    Vector var1(rows, 0.0), var2(rows, 0.0);
    Vector se1(rows, 0.0), se2(rows, 0.0);

    for (size_t i = 0; i < rows; ++i) {
        Vector values1, values2, weights;
        for (const auto& matrix : matrices) {
            values1.push_back(matrix[i][0]);
            values2.push_back(matrix[i][1]);
            weights.push_back(matrix[i][2]);
        }

        // Compute weighted mean
        mu1[i] = calculate_weighted_mean(values1, weights);
        mu2[i] = calculate_weighted_mean(values2, weights);

        // Compute weighted variance
        var1[i] = calculate_weighted_variance(values1, weights, mu1[i]);
        var2[i] = calculate_weighted_variance(values2, weights, mu2[i]);

        // Compute weighted SEM
        se1[i] = calculate_weighted_sem(var1[i], weights);
        se2[i] = calculate_weighted_sem(var2[i], weights);
    }

    // Save mu1, mu2, var1, var2, se1, se2 to a file
    // 6 columns for mu1, mu2, var1, var2, se1, se2
    Matrix final(rows, Vector(6, 0.0));  
    for (size_t i = 0; i < rows; ++i) {
        final[i][0] = mu1[i];
        final[i][1] = mu2[i];
        final[i][2] = var1[i];
        final[i][3] = var2[i];
        final[i][4] = se1[i];
        final[i][5] = se2[i];
    }

    save_matrices_to_file(save_path, final);
}

void cache_size(const std::vector<std::string>& all_filenames, const std::string& substring, const std::string& save_path) {
    std::vector<Matrix> matrices;
    for (const auto& filename : all_filenames) {
        if (filename.find(substring) != std::string::npos) {
            Matrix mat = read_file(filename);
            if (mat.empty() || mat[0].empty()) {
                std::cerr << "Warning: Empty or invalid matrix read from file: " << filename << std::endl;
                continue;
            }
            matrices.push_back(mat);
        }
    }
    
    if (matrices.empty()) {
        std::cerr << "Error: No valid matrices found for substring: " << substring << std::endl;
        return;
    }
    
    size_t rows = matrices[0].size();
    size_t cols = matrices[0][0].size();

    // Check for consistent matrix dimensions
    for (const auto& matrix : matrices) {
        if (matrix.size() != rows || matrix[0].size() != cols) {
            std::cerr << "Error: Inconsistent matrix dimensions." << std::endl;
            return;
        }
    }

    // Normalize by the first element of each row
    for (auto& matrix : matrices) {
        for (size_t i = 0; i < rows; ++i) {
            double n = matrix[i][0];
            if (n == 0) {
                std::cerr << "Error: Division by zero detected during normalization for row " << i << std::endl;
                return;
            }
            for (size_t j = 0; j < cols; ++j) {
                matrix[i][j] /= n;
            }
        }
    }

    // Calculate the mean, standard deviation, and SEM
    Matrix mu(rows, Vector(cols, 0.0));
    Matrix sd(rows, Vector(cols, 0.0));
    Matrix stderr(rows, Vector(cols, 0.0));

    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            Vector column_values;
            for (const auto& matrix : matrices) {
                column_values.push_back(matrix[i][j]);
            }

            mu[i][j] = calculate_mean(column_values);

            sd[i][j] = calculate_std_dev(column_values, mu[i][j]);

            stderr[i][j] = calculate_sem(sd[i][j], matrices.size());
        }
    }

    save_matrices_to_file(save_path, mu);
}

int main() {
    std::string RESULTS_DIRECTORY = "data";
    std::string FINAL_DIRECTORY = "final";

    // Check if RESULTS_DIRECTORY exists and is a directory
    if (!fs::exists(RESULTS_DIRECTORY) || !fs::is_directory(RESULTS_DIRECTORY)) {
        std::cerr << "Error: Directory " << RESULTS_DIRECTORY << " does not exist or is not a directory." << std::endl;
        return 1; 
    }

    // Attempt to create FINAL_DIRECTORY; check for failure
    std::error_code ec;
    if (!fs::create_directory(FINAL_DIRECTORY, ec) && ec) {
        std::cerr << "Error: Could not create directory " << FINAL_DIRECTORY << ". Error code: " << ec.message() << std::endl;
        return 1;
    }

    std::vector<std::string> filenames = get_all_results_filenames(RESULTS_DIRECTORY);

    // Check if any filenames were obtained
    if (filenames.empty()) {
        std::cerr << "Error: No files found in directory " << RESULTS_DIRECTORY << std::endl;
        return 1;
    }

    obs1(filenames, "_energy.txt", FINAL_DIRECTORY + "/energy.txt");
    obs1(filenames, "_energy_IS.txt", FINAL_DIRECTORY + "/energy_IS.txt");
    ridge(filenames, "_ridge_E.txt", FINAL_DIRECTORY + "/ridge_E.txt");
    ridge(filenames, "_ridge_S.txt", FINAL_DIRECTORY + "/ridge_S.txt");
    obs1(filenames, "_acceptance_rate.txt", FINAL_DIRECTORY + "/acceptance_rate.txt");
    obs1(filenames, "_inherent_structure_timings.txt", FINAL_DIRECTORY + "/inherent_structure_timings.txt");
    obs1(filenames, "_walltime_per_waitingtime.txt", FINAL_DIRECTORY + "/walltime_per_waitingtime.txt");
    cache_size(filenames, "_cache_size.txt", FINAL_DIRECTORY + "/cache_size.txt");

    return 0; 
}
