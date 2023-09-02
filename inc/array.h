#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#ifndef HDSPIN_ARRAY_H
#define HDSPIN_ARRAY_H


namespace hdspin_array
{

template<typename T>
class Dynamic2DArray
{
private:
    std::vector<std::vector<T>> data;
    int rows;
    int cols;

public:
    Dynamic2DArray(){};

    std::vector<int> shape() const
    {
        return std::vector<int>{rows, cols};
    }

    void load(const std::string& filename, const std::string& delimiter)
    {
        std::ifstream file(filename);

        if (!file.is_open())
        {
            throw std::runtime_error("Unable to open file: " + filename);
        }

        data.clear();
        std::string line;
        while (std::getline(file, line))
        {
            std::vector<T> row;
            std::istringstream iss(line);
            std::string value;
            while (std::getline(iss, value, delimiter.c_str()[0]))
            {
                row.push_back(std::stoi(value));
            }
            data.push_back(row);
        }
        file.close();

        if (data.empty()) {
            throw std::runtime_error("No data found in file.");
        }

        cols = data[0].size();
        rows = data.size();
    }

    void print() const
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }


};


}

#endif
