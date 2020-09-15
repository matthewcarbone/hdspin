executables := $(shell find . -maxdepth 1 -name "*.out")

# rr:
# 	mpic++ -fopenmp hdspin/main.cpp -std=c++11 -o main.out hdspin/gillespie.cpp hdspin/standard.cpp hdspin/utils/general_utils.cpp hdspin/utils/init_utils.cpp

# Run with, e.g., mpirun -np 1 ./main.out test/file.txt 100 8 0.75 1.0 0 1 100
# local:
# 	/usr/local/bin/mpic++ -Xpreprocessor -fopenmp -lomp hdspin/main.cpp -std=c++11 -o main.out hdspin/gillespie.cpp hdspin/standard.cpp hdspin/utils/general_utils.cpp hdspin/utils/init_utils.cpp

rr: 
	g++ -Wall hdspin/main.cpp -fopenmp -std=c++11 -o main.out hdspin/gillespie.cpp hdspin/standard.cpp hdspin/utils/general_utils.cpp hdspin/utils/init_utils.cpp hdspin/utils/grid_utils.cpp

local:
	g++ -Wall hdspin/main.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -o main.out hdspin/gillespie.cpp hdspin/standard.cpp hdspin/utils/general_utils.cpp hdspin/utils/init_utils.cpp hdspin/utils/grid_utils.cpp

clean:
	-rm $(executables)
