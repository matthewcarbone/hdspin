executables := $(shell find exe -maxdepth 1 -name "*.out")

INC=-I inc
SRC=$(wildcard src/*/*.cpp)

CC=mpic++

# Previous local flags using -fopenmp
# LOCAL_FLAGS=-Xpreprocessor -fopenmp -lomp -std=c++17
LOCAL_FLAGS=-std=c++17

linux:
	$(CC) $(INC) $(LOCAL_FLAGS) main/main.cpp -o exe/main.out $(SRC) -O3 

debug:
	$(CC) $(INC) $(LOCAL_FLAGS) main/debug.cpp -o exe/debug.out $(SRC)

test:
	$(CC) $(INC) $(LOCAL_FLAGS) main/test.cpp -o exe/test.out $(SRC)

clean:
	-rm $(executables)
