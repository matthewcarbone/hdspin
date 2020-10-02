#ifndef OBS_BASE_H
#define OBS_BASE_H

#include <iostream>
#include <vector>
#include <fstream>

class Base
{
protected:
    FILE *outfile;
    int grid_length;

public:
    void open_outfile(const std::string);
    void close_outfile();  // Close the output stream when done
};

#endif
