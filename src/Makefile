#Generic c++ Makefile

CXX = g++
CXXFLAGS = -O3 -g -fopenmp
LINKS = 

MAIN = cntGenerate

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

$(MAIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LINKS)

.PHONY: clean
clean:
	rm -f $(OBJ) $(MAIN)
