executables := $(shell find exe -maxdepth 1 -name "*.out")

INC=-I inc
SRC = $(wildcard src/*/*.cpp)

CC=mpic++
LOCAL_CC=/usr/local/bin/mpic++
LINUX_LOCAL_CC=/usr/bin/mpic++

# Previous local flags using -fopenmp
# LOCAL_FLAGS=-Xpreprocessor -fopenmp -lomp -std=c++17
LOCAL_FLAGS=-std=c++17

# Previous remote flags
# REMOTE_FLAGS_1=-fopenmp -std=c++17
REMOTE_FLAGS_1=-std=c++17

rr: 
	$(CC) $(INC) $(REMOTE_FLAGS_1) main/main.cpp -o exe/main.out $(SRC) -O3

local: 
	$(LOCAL_CC) $(INC) $(LOCAL_FLAGS) main/main.cpp -o exe/main.out $(SRC) -O3

local_linux:
	$(LINUX_LOCAL_CC) $(INC) $(LOCAL_FLAGS) main/main.cpp -o exe/main.out $(SRC) -O3 

debug:
	$(CC) $(INC) $(LOCAL_FLAGS) main/debug.cpp -o exe/debug.out $(SRC)

test:
	$(CC) $(INC) $(LOCAL_FLAGS) main/test.cpp -o exe/test.out $(SRC)

clean:
	-rm $(executables)
