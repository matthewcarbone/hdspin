#ifndef OBS_BASE_H
#define OBS_BASE_H

#include <iostream>
#include <vector>
#include <fstream>

#include "Utils/structures.h"

class Base
{
protected:
    FileNames fnames;
    FILE *outfile;
public:
    Base(const FileNames);
    ~Base();
};

#endif
