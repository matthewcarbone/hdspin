executables := $(shell find . -maxdepth 1 -name "*.out")

rr:
	module load compilers/gcc-10.1.0
	mpic++ -fopenmp hdspin/main.cpp -std=c++11 -o main.out hdspin/utils/general_utils.cpp

# Run with, e.g., mpirun -np 1 ./main.out test/file.txt 100 8 0.75 1.0 0 1 100
local:
	/usr/local/bin/mpic++ -Xpreprocessor -fopenmp -lomp hdspin/main.cpp -std=c++11 -o main.out hdspin/utils/general_utils.cpp

clean:
	-rm $(executables)
