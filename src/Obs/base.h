#include "Obs/base.h"

void Base::open_outfile(const std::string d)
{
    outfile = fopen(d.c_str(), "w");
}

void Base::close_outfile(){fclose(outfile);}
