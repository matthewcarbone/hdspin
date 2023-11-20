// ChatGPT 3.5 was used to assist in the creation of this cpp file

// Could I have used an existing array library to do this? Absolutely. But
// I wanted to feel the pain this time, and learn how to work with C++'s
// "beautiful" filesystem interface.

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

std::vector<std::string> get_matching_files(const std::string& directoryPath, const std::string& pattern) {
    std::vector<std::string> matchingFiles;

    try {
        for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
            bool is_regular_path = std::filesystem::is_regular_file(entry.path());
            bool npos_condition = std::filesystem::path(entry).filename().string().find(pattern) != std::string::npos;
            if (is_regular_path && npos_condition) {
                matchingFiles.push_back(entry.path().string());
            }
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error accessing directory: " << ex.what() << std::endl;
    }

    std::sort(matchingFiles.begin(), matchingFiles.end());

    return matchingFiles;
}


// int main() {
//     // Replace "your_directory_path" with the actual directory path
//     std::string directoryPath = "your_directory_path";
//     // Replace "_test.txt" with the actual pattern
//     std::string pattern = "_test.txt";

//     std::vector<std::string> matchingFiles = getMatchingFiles(directoryPath, pattern);

//     if (!matchingFiles.empty()) {
//         std::cout << "Matching Files:" << std::endl;
//         for (const auto& file : matchingFiles) {
//             std::cout << file << std::endl;
//         }
//     } else {
//         std::cout << "No matching files found." << std::endl;
//     }

//     return 0;
// }



template <typename T>
class Array2D {
private:
    std::vector<std::vector<T>> data;

public:
    void load(const std::string fileName)
    {
        std::ifstream inFile(fileName);
        if (!inFile) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            exit(1);
        }

        std::string line;
        while (std::getline(inFile, line)) {
            std::istringstream iss(line);
            T value;
            data.emplace_back(); // Add a new row

            while (iss >> value) {
                data.back().push_back(value);
            }
        }

        inFile.close();
    }

    void display() const {
        for (const auto& row : data) {
            for (const auto& element : row) {
                std::cout << element << " ";
            }
            std::cout << std::endl;
        }
    }

    size_t getNumRows() const {
        return data.size();
    }

    size_t getNumCols() const {
        return (data.empty() ? 0 : data[0].size());
    }
};


int main() {
    // Replace "your_filename.txt" with the actual file name
    // std::string fileName = "test.txt";

    // // Example with integer type
    // Array2D<double> arr;
    // arr.load(fileName);

    // std::cout << "Integer Array:" << std::endl;
    // arr.display();

    // std::cout << "Number of Rows: " << arr.getNumRows() << std::endl;
    // std::cout << "Number of Columns: " << arr.getNumCols() << std::endl;

    // You can repeat the process for other data types as needed.

    const std::string directoryPath = "data/";
    const std::string pattern = "_energy.txt";

    std::vector<std::string> matchingFiles = get_matching_files(directoryPath, pattern);

    // if (!matchingFiles.empty()) {
    //     std::cout << "Matching Files:" << std::endl;
    //     for (const auto& file : matchingFiles) {
    //         std::cout << file << std::endl;
    //     }
    // } else {
    //     std::cout << "No matching files found." << std::endl;
    // }

    std::vector<Array2D<float>> loaded_arrays;
    for (const auto& file : matchingFiles) {
        Array2D<float> arr;
        arr.load(file);
        loaded_arrays.push_back(arr);
    }

    return 0;
}
