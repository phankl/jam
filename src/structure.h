#ifndef STRUCTURE_HEADER
#define STRUCTURE_HEADER

#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <numbers>

#include <complex>
#include <fftw3.h>

#include "analysis.h"
#include "tube.h"
#include "xyz.h"

using namespace std;

class Structure{
  
  public:
   
    Structure(int, int, int, double, double, XYZ, XYZ, XYZ);
    
    vector<double> odf(int);
    void printDataFile(string);

  private:

    int nTube;
    double rTube;
    double lTube;

    XYZ boxSize;

    XYZ mean;
    XYZ std;

    vector<Tube> tubes;
    vector<Tube> ghostTubes;

    vector<XYZ> atoms;
    
    double distance(Tube, Tube, double, double);
    bool checkOverlap(Tube);

    void generate(int);
    void segment(int);
};

#endif
