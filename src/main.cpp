#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>

#include "structure.h"
#include "xyz.h"

using namespace std;

int main() {
 
  // structure properties

  double boxLength = 500.0;
  double volumeFraction = 0.05;

  int seed = 2022;
  int nSegment = 10;
  XYZ boxSize(boxLength, boxLength, boxLength);
  XYZ mean(0.0, 0.0, 5.0);
  XYZ std(1, 1, 1);

  // CNT properties
  
  double aCC = 1.421e-1;
  double mC = 12.0;

  int n = 10;
  int m = 10;
  double lTube = 100.0;
  
  double rTube = 0.5 * sqrt(3.0*(n*n + n*m + m*m)) / numbers::pi * aCC;
  int nTube = volumeFraction * pow(boxLength, 3) / (numbers::pi * pow(rTube, 2) * lTube);
  
  double linearDensity = 2.0 / 3.0 * sqrt(n*n + n*m + m*m) * mC / aCC;
  double mTube = linearDensity * lTube;

  double cutoff = 1.1 * lTube / nSegment;
  
  cout << "Generating structure with " << nTube << " CNTs at volume fraction of " << 100*volumeFraction << "%." << endl;

  Structure structure(seed, nTube, nSegment, rTube, lTube, mTube, boxSize, mean, std);
  structure.printDataFile("cnt.data");
  structure.printInputFile(cutoff, cutoff, "cnt.in");

  return 0;
}
