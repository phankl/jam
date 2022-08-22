#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "structure_rigid.h"
#include "structure_elastic.h"
#include "structure_film.h"
#include "xyz.h"

using namespace std;

int main() {

  double pi = 4 * atan(1);
 
  // structure properties
    
  int seed = 2023;
  XYZ boxSize(20000.0, 20000.0, 200.0);
  XYZ mean(0.0, 0.0, 0.0);
  XYZ std(1, 1, 10);
  
  double density = 0.14 / 1.6605;

  /*
    alignment paper values
    int seed = 2022;
    XYZ boxSize(10000.0, 10000.0, 10000.0);
    XYZ mean(0.0, 0.0, 1.0);
    XYZ std(1, 1, 3);

    double density = 0.002 / 1.6605;
  */

  // CNT properties
  
  double aCC = 1.421;
  double mC = 12.0;

  int n = 10;
  int m = 10;
  int nSegment = 500;
  double lTube = 10000.0;
  
  double rTube = 0.5 * sqrt(3.0*(n*n + n*m + m*m)) / pi * aCC;
  double vBox = boxSize[0] * boxSize[1] * boxSize[2];
  
  double linearDensity = 4.0 / 3.0 * sqrt(n*n + n*m + m*m) * mC / aCC;
  double mTube = linearDensity * lTube;
  
  int nTube = density * vBox / mTube;

  double skin = 2 * aCC;
 
  // simulation properties
  
  int steps = 1000000;
  double temp = 600.0;

  cout << "Mass density: " << 1.6605 * density << endl;

  // StructureFilm structure(seed, nTube, nSegment, rTube, lTube, mTube, skin, boxSize, mean, std);
  StructureFilm structure(seed, nTube, nSegment, rTube, lTube, mTube, skin, boxSize, mean, std);
  structure.printDataFile("data.iso");
  // structure.printInertiaFile("cnt.inertia");
  // structure.printInputFile(steps, temp, "cnt.in");

  int nBins = 50;
  vector<double> odf2D = structure.odf2D(nBins, {1, 0, 0}, {0, 0, 1});
  vector<double> odf3D = structure.odf3D(nBins, {0, 0, 1});

  ofstream odf2DFile("odf2D.dat");
  ofstream odf3DFile("odf3D.dat");

  for (int i = 0; i < nBins + 1; i++) {
    double x = i * pi / nBins;
    if (i == nBins) {
      odf2DFile << x << " " << odf2D[0] << endl;
      odf3DFile << x << " " << odf3D[0] << endl;
    } 
    else {
      odf2DFile << x << " " << odf2D[i] << endl;
      odf3DFile << x << " " << odf3D[i] << endl;
    }
  }

  odf2DFile.close();
  odf3DFile.close();

  return 0;
}
