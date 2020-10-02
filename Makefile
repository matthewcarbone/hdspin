executables := $(shell find exe -maxdepth 1 -name "*.out")

INC=-I inc
SRC = $(wildcard src/*/*.cpp)

CC=g++
LOCAL_FLAGS=-Xpreprocessor -fopenmp -lomp -std=c++17
REMOTE_FLAGS_1=-fopenmp -std=c++17

# rr:
# 	mpic++ -fopenmp hdspin/main.cpp -std=c++11 -o main.out hdspin/gillespie.cpp hdspin/standard.cpp hdspin/utils/general_utils.cpp hdspin/utils/init_utils.cpp

# Run with, e.g., mpirun -np 1 ./main.out test/file.txt 100 8 0.75 1.0 0 1 100
# local:
# 	/usr/local/bin/mpic++ -Xpreprocessor -fopenmp -lomp hdspin/main.cpp -std=c++11 -o main.out hdspin/gillespie.cpp hdspin/standard.cpp hdspin/utils/general_utils.cpp hdspin/utils/init_utils.cpp

rr: 
	$(CC) $(INC) $(REMOTE_FLAGS_1) main/main.cpp -o exe/main.out $(SRC) -O3

debug:
	$(CC) $(INC) $(LOCAL_FLAGS) main/debug.cpp -o exe/debug.out $(SRC)

test:
	$(CC) $(INC) $(LOCAL_FLAGS) main/test.cpp -o exe/test.out $(SRC)

clean:
	-rm $(executables)
