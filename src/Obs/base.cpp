#include "Obs/base.h"
#include "Utils/structures.h"

Base::Base(const FileNames fnames) : fnames(fnames) {};

Base::~Base(){fclose(outfile);}
