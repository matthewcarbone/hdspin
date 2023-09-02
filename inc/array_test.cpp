#include "array.h"

int main() {
    hdspin_array::Dynamic2DArray<double> arr;

    std::string filename = "test.txt";
    arr.load(filename, " ");
    arr.print();


    // try {
        
    //     double rowMean = arr.mean(0);
    //     double colMean = arr.mean(1);
    //     double rowStdDev = arr.standardDeviation(0);
    //     double colStdDev = arr.standardDeviation(1);

    //     std::cout << "Row Mean: " << rowMean << std::endl;
    //     std::cout << "Column Mean: " << colMean << std::endl;
    //     std::cout << "Row Standard Deviation: " << rowStdDev << std::endl;
    //     std::cout << "Column Standard Deviation: " << colStdDev << std::endl;

    // } catch (const std::exception& e) {
    //     std::cerr << "Error: " << e.what() << std::endl;
    // }

    return 0;
}
