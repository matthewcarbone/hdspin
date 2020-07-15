executables := $(shell find . -maxdepth 1 -name "*.o")

#export PATH=/opt/pb/gcc-10.1.0/bin:$PATH
#export LD_LIBRARY_PATH=/opt/pb/gcc-10.1.0/lib:/opt/pb/gcc-10.1.0 lib64:$LD_LIBRARY_PATH

all:
	g++ -std=c++17 -Wall main.cpp -o main.o hdspin/gillespie.cpp \
	hdspin/utils/grid_utils.cpp hdspin/utils/init_utils.cpp \
	hdspin/utils/general_utils.cpp hdspin/utils/file_utils.cpp

clean:
	-rm $(executables)

