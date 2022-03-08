#Generic c++ Makefile

CXX = g++-10
CXXFLAGS = -O3 -g -fopenmp -std=c++2a
LINKS = 

MAIN = cntGenerate

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

$(MAIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LINKS)

.PHONY: clean
clean:
	rm -f $(OBJ) $(MAIN)
