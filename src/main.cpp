#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>

#include "structure.h"
#include "xyz.h"

using namespace std;

int main() {
 
  int seed = 2022;
  int nTube = 1000;
  int nSegment = 10;
  double rTube = 0.6785;
  double lTube = 1000.0;
  XYZ boxSize(10000.0, 10000.0, 10000.0);
  XYZ mean(0.0, 0.0, 0.0);
  XYZ std(1, 1, 1);

  Structure structure(seed, nTube, nSegment, rTube, lTube, boxSize, mean, std);
  structure.printDataFile("cnt.data");
  
  return 0;
}
