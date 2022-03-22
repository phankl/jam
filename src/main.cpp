#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "structure.h"
#include "xyz.h"

using namespace std;

int main() {

  double pi = 4 * atan(1);
 
  // structure properties

  double boxLength = 5000.0;
  double volumeFraction = 0.01;

  int seed = 2021;
  int nSegment = 50;
  XYZ boxSize(boxLength, boxLength, boxLength);
  XYZ mean(0.0, 0.0, 5.0);
  XYZ std(1, 1, 1);

  // CNT properties
  
  double aCC = 1.421;
  double mC = 12.0;

  int n = 10;
  int m = 10;
  double lTube = 1000.0;
  
  double skin = 1.0;
  double rTube = 0.5 * sqrt(3.0*(n*n + n*m + m*m)) / pi * aCC;
  int nTube = volumeFraction * pow(boxLength, 3) / (pi * pow(rTube, 2) * lTube);
  
  double linearDensity = 2.0 / 3.0 * sqrt(n*n + n*m + m*m) * mC / aCC;
  double mTube = linearDensity * lTube;

  double cutoff = 3.5 * lTube / nSegment;
  
  cout << "Generating structure with " << nTube << " CNTs at volume fraction of " << 100*volumeFraction << "%." << endl;

  Structure structure(seed, nTube, nSegment, rTube+skin, lTube, mTube, boxSize, mean, std);
  structure.printDataFile("cnt.data");
  structure.printInertiaFile("cnt.inertia");
  // structure.printInputFile(cutoff, cutoff, "cnt.in");

  return 0;
}
